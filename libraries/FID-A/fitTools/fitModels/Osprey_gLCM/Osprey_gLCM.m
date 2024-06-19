function [ModelParameter] = Osprey_gLCM(DataToModel, JsonModelFile, average, noZF, scaleData, NumericJacobian, CheckGradient, BasisSetStruct)
%% Global function for new Osprey LCM
% Inputs:   DataToModel     - FID-A/Osprey struct with data or cell of structs
%           JsonModelFile   - Master model file for all steps
%           average         - average data prior to modeling (NIfTI-MRS only)
%           noZF            - avoid zero-filling
%           scaleData       - Do data to basis scaling or apply scale factor
%           NumericJacobian - use numerical jacobian flag
%           CheckGradient   - Do a gradient check in the lsqnonlin solver
%           BasisSetStruct  - Include predefined basisset struct
% Outputs:  explicit        - Struct with model parameters
%           implicit        - NII results file

%% 0. Check inputs
arguments
    % Argument validation introduced in relatively recent MATLAB versions
    % (2019f?)
    DataToModel {isStructOrCell}
    JsonModelFile {isCharOrStruct}
    average double {mustBeNumeric} = 0;             % optional
    noZF double {mustBeNumeric} = 0;                % optional
    scaleData double {mustBeNumeric} = 0;           % optional
    NumericJacobian double {mustBeNumeric} = 0;     % optional
    CheckGradient double {mustBeNumeric} = 0;       % optional
    BasisSetStruct struct = [];                     % optional
end

% If input is just a single struct, move it into a cell array
if ~iscell(DataToModel)                                                       % Just a single file
    temporaryCell       = {};
    temporaryCell{1}    = DataToModel;                                        % Default input is a cell array so we need to change it
    DataToModel         = temporaryCell;
else
    % If input is a cell, check whether it is a NII file name that can be
    % directly parsed
    if isstring(DataToModel{1})
        if contains(DataToModel{1},'.nii')                                    % We do only accept .nii files for direct parsing
            DataToModel{1}             = io_loadspec_niimrs(DataToModel{1});  % Load data
        end
    end
end

%% What happens here:
%   1. Decode model json file
%   2. Prepare basisset matrix (and export as NII)
%       2a. Do on-the-fly generation of MMs
%   3. Prepare data according to model json & Run steps defined in model json
%   4. Save results (and export as NII)

%% 1. Decode model json file

% Read the json file and generate a ModelProcedure struct from it which will guide the
% rest of the analysis. Catch missing parameters here?
if ~isstruct(JsonModelFile)                                                   %Parsed as json file
    ModelProcedure = jsonToStruct(JsonModelFile);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
else
    ModelProcedure  = JsonModelFile;
end

% Check if the R part is used for optimization and a factor 2 zerofilling
zf = 0;
for ss = 1 : length(ModelProcedure.Steps)
    if ~isfield(ModelProcedure.Steps{ss},'fit_opts')
        ModelProcedure.Steps{ss}.fit_opts.optimSignalPart = 'R';
    end
    if isfield(ModelProcedure.Steps{ss}.fit_opts,'optimSignalPart')
        if strcmp(ModelProcedure.Steps{ss}.fit_opts.optimSignalPart,'R')
            zf = 1;
        end
    end
end
if noZF                                                                     % Overwrite zero-fill
    zf = 0;
end

%% 2. Prepare basisset matrix (and export as NII)
% Load basisset files, add MMs if needed, resample basis sets according to
% the DataToModel. Generate a basisset matrix for each step? including the
% indirect dimensions for MSM.
if isempty(BasisSetStruct)                                                  % User supplied a recalculated basis set
    if length(ModelProcedure.basisset.file) == 1
        basisSet = load(ConvertRelativePath(ModelProcedure.basisset.file{1})); % Load basis set
        basisSet = basisSet.BASIS;
        basisSet = recalculateBasisSpecs(basisSet);                         % Add ppm axis and frequency domain data
        basisSet = fit_sortBasisSet(basisSet);                              % Sort according to Osprey standard
        if isfield(ModelProcedure.basisset,'opts')                          % Apply options to basis, e.g. pick subspectra
            if isfield(ModelProcedure.basisset.opts,'index')               % Index for subspectra
                basisSet.fids = basisSet.fids(:,:,ModelProcedure.basisset.opts.index);  % Pick fids according to index
                basisSet.specs = basisSet.specs(:,:,ModelProcedure.basisset.opts.index); % Pick specs according to index
                basisSet.sz = size(basisSet.fids);                              % Recalculate size entry
                try
                    basisSet.nExtra = basisSet.sz(3);                               % Update extra dimension
                catch
                end
            end
        end
    else
        for bb = 1 : length(ModelProcedure.basisset.file)
            if bb == 1
                basisSet = load(ModelProcedure.basisset.file{bb});          % Load basis set
                basisSet = basisSet.BASIS;
                basisSet = recalculateBasisSpecs(basisSet);                 % Add ppm axis and frequency domain data
                basisSet = fit_sortBasisSet(basisSet);                      % Sort according to Osprey standard
            else
                basisSetToAdd = load(ModelProcedure.basisset.file{bb});     % Load basis set
                basisSetToAdd = basisSetToAdd.BASIS;
                basisSetToAdd = recalculateBasisSpecs(basisSetToAdd);       % Add ppm axis and frequency domain data
                basisSetToAdd = fit_sortBasisSet(basisSetToAdd);            % Sort according to Osprey standard
                basisSet.fids = cat(3,basisSet.fids,basisSetToAdd.fids);    % Concatenate time domain basis functions for 2D fit
                basisSet.specs = cat(3,basisSet.specs,basisSetToAdd.specs); % Concatenate frequency domain basis functions for 2D fit
            end
        end
        basisSet.sz = size(basisSet.fids);                                  % Recalculate size entry
        basisSet.nExtra = basisSet.sz(3);                                   % Update extra dimension
    end
    if isfield(ModelProcedure.basisset,'opts')                              % Apply options to basis, e.g. repeat for averages
        if isfield(ModelProcedure.basisset.opts,'repeat')                   % Repeat for each average
            basisSet.fids = repmat(basisSet.fids,[1 1 ModelProcedure.basisset.opts.repeat]);    % Repeat time domain basis functions for 2D fit
            basisSet.specs = repmat(basisSet.specs,[1 1 ModelProcedure.basisset.opts.repeat]);  % Repeat frequency domain basis functions for 2D fit
            basisSet.sz = size(basisSet.fids);                              % Recalculate size entry
            basisSet.nExtra = basisSet.sz(3);                               % Update extra dimension
        end
    end

else
    basisSet = BasisSetStruct;                                              % Take user supplied basis set
    basisSet.nExtra = 0;                                                    % Normally has no extra dimension
end

if length(scaleData) > 1                                                    % The user has supplied a scaling factor for the data               
    DoDataScaling = 2;                                                      % Apply scale to data
else if scaleData ~= 0                                                      % Do we want to perform data scaling?
        DoDataScaling = 1;                                                  % Apply scale to data
        if scaleData ~= 1
            DoDataScaling = 2;
        end
    else
        DoDataScaling = 0;                                                  % No data scaling
    end
end

%% 2.a Generate MM/Lip basis functions
if isfield(ModelProcedure.basisset, 'mmdef') && ~isempty(ModelProcedure.basisset.mmdef)                    %Overwrite if defined
    ModelProcedure.basisset.mmdef{1} = 'which(fullfile(''libraries'',''FID-A'',''fitTools'',''fitModels'',''Osprey_gLCM'',''fitClass'',''mm-definitions'',''MMLipLCModel.json''))';

    % First, we'll remove any MM or lipid basis functions that may be in the
    % input basis set (to avoid duplication)
    basisSet = scrubMMLipFromBasis(basisSet);
    
    % Read out MM/Lip configuration file and create matching
    MMLipConfig = jsonToStruct(ConvertRelativePath(ModelProcedure.basisset.mmdef{1}));
    [basisSim] = makeMMLipBasis(basisSet, MMLipConfig, DataToModel{1});
    % Join basis sets together
    basisSet = joinBasisSets(basisSet, basisSim);
end

%% 3. Prepare data according to model json and model data
% Prepare data for each fit step, again, including the indirect dimensions
% for MSM
% Create the spline basis functions for the given resolution, fit range,
% and knot spacing parameter.

% Cell array of data which we will loop over
for kk = 1 : length(DataToModel)
    if ~isstruct(DataToModel{kk})                                                  % It is not a FID-A struct?
        if contains(DataToModel{kk},'.nii')                                        % We do only accept .nii files for direct parsing
            temp                = io_loadspec_niimrs(DataToModel{kk});          % Temporarily load the NIfTI-MRS file
            if average                                                        % Average before modeling
                temp            = op_averaging(temp);                         % Average before modeling
            else
                temp.dims.extras = temp.dims.averages;
                temp.dims.averages = 0;                                       % Multi-spectral modelign only works in extra dim
            end
            DataToModel{kk} 	= temp;
        end
    end

    if DoDataScaling == 2
        DataToModel{kk}   = op_ampScale(DataToModel{kk}, 1/scaleData(kk));                     % apply scale to data 
    else if DoDataScaling == 1
         scaleData(kk) = max(real(DataToModel{kk}.specs(DataToModel{kk}.ppm > -2 & DataToModel{kk}.ppm < 10 ,:))) / max(max(max(real(basisSet.specs(basisSet.ppm > -2 & basisSet.ppm < 10 ,:)))));
         DataToModel{kk}   = op_ampScale(DataToModel{kk}, 1/scaleData(kk));                    % apply scale to data
    end
    end

    % Apply zero-filling if needed
    if zf &&  ~DataToModel{kk}.flags.zeropadded
        DataToModel{kk} = op_zeropad(DataToModel{kk},2);                    % Zero-fill data if needed (real part optimization)
    end

    for ss = 1 : length(ModelProcedure.Steps)
        opts = [];                                                                      % we have to clean older steps out
        if isfield(ModelProcedure.Steps{ss},'fit_opts')
            opts  = ModelProcedure.Steps{ss}.fit_opts;
        end
        % Setup options from model procedure json file
        if isfield(ModelProcedure.Steps{ss},'ModelFunction')
            opts.ModelFunction      = ModelProcedure.Steps{ss}.ModelFunction;           % get model function
        else
            opts.ModelFunction      = 'GeneralizedBasicPhysicsModel';                   % set default model function
        end
        if isfield(ModelProcedure.Steps{ss},'parameter')
            opts.parameter  = ModelProcedure.Steps{ss}.parameter;                       % specify parametrizations constructor
        else
            opts.parameter  = {'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'}; % set default parametrizations constructor
        end

        if ~isfield(opts,'solver')
            opts.solver             = 'lsqnonlin';                                  % set default solver
        end       
        if isfield(ModelProcedure.Steps{ss},'parametrizations')
            opts.parametrizations  = ModelProcedure.Steps{ss}.parametrizations;         % specify parametrizations constructor
            parameter = {'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'};
            opts.parametrizations  = orderfields(opts.parametrizations,parameter(ismember(parameter, fieldnames(opts.parametrizations)))); % order struct names according to standard
        end
        if ~isfield(ModelProcedure.Steps{ss},'basisset')                       % which basis functions to include
            ModelProcedure.Steps{ss}.basisset.spec = 1;
            ModelProcedure.Steps{ss}.basisset.include = {'Asc','Asp','Cr','CrCH2','GABA','GPC','GSH','Gln','Glu','mI','Lac','NAA',...
                                                        'NAAG','PCh','PCr','PE','sI','Tau','Lip13a','Lip13b','Lip13c','Lip13d', 'Lip09','MM09','Lip20','MM12','MM14','MM17','MM20'};
        else
            if ~isfield(ModelProcedure.Steps{ss}.basisset, 'include')
            ModelProcedure.Steps{ss}.basisset.include = {'Asc','Asp','Cr','CrCH2','GABA','GPC','GSH','Gln','Glu','mI','Lac','NAA',...
                                                         'NAAG','PCh','PCr','PE','sI','Tau','Lip13a','Lip13b','Lip13c','Lip13d', 'Lip09','MM09','Lip20','MM12','MM14','MM17','MM20'};
            end
            if ~isfield(ModelProcedure.Steps{ss}.basisset, 'spec')
                ModelProcedure.Steps{ss}.basisset.spec = 1;
            end
        end
        if isfield(DataToModel{kk},'FWHM') && ss == 1                                   % get inital FWHM estimate 
            opts.parametrizations.gaussLB.init = DataToModel{kk}.FWHM;
        end
        if isfield(ModelProcedure.Steps{ss}, 'extra')                                   % Is 2D model
            if ModelProcedure.Steps{ss}.extra.flag == 1 && isfield(ModelProcedure.Steps{ss}.extra,'DynamicModelJson')
                if ~isstruct(ModelProcedure.Steps{ss}.extra.DynamicModelJson)                                               %Parsed as json file
                    opts.paraIndirect  = jsonToStruct(ConvertRelativePath(ModelProcedure.Steps{ss}.extra.DynamicModelJson)); % specify parametrizations constructor for indirect dimension
                else
                    opts.paraIndirect = ModelProcedure.Steps{ss}.extra.DynamicModelJson;
                end
            end
        end
        if strcmp(ModelProcedure.Steps{ss}.module,'OptimReg')                           % Change tolerance values (good for regularization)
            opts.FunctionTolerance = 1e-3;
            opts.StepTolerance = 1e-3;
            opts.OptimalityTolerance = 1e-3;
        else
            opts.FunctionTolerance = 1e-10;
            opts.StepTolerance = 1e-10;
            opts.OptimalityTolerance = 1e-10;
        end
        opts.NumericJacobian    = NumericJacobian;                                      % Use numerical jacobian instead of the analytical jacobian
        opts.CheckGradient      = CheckGradient;                                        % Do a gradient check in the lsqnonlin solver

        clc
        fprintf('Running model procedure step %i. \n', ss);

        % Create an instance of the class
        if ss == 1
            ModelParameter{kk,1} = FitObject(DataToModel{kk}, basisSet, opts);  % Create OspreyFitObj instance (ss = 1)
        else
            ModelParameter{kk,1}.updateOptsAccordingToStep(opts);               % Update model options according to step
        end

        % Set basisset according to step
        ModelParameter{kk,1}.excludeBasisFunctionFromFit('all');
        ModelParameter{kk,1}.includeBasisFunctionInFit(ModelProcedure.Steps{ss}.basisset.include);

        % Set the first set of expectation and standard deviation values
        % for the non-linear parameters
        if isfield(ModelProcedure.basisset, 'mmdef') && ~isempty(ModelProcedure.basisset.mmdef)
            % Add the initialization of EX/SD values (for frequency shift and
            % Lorentzian LB) here:
            [EXT2, SDT2, SDSH] = fit_setExSDValues(DataToModel{kk}, basisSet, MMLipConfig, 0);
             % Store EXT2, SDT2, SDSH to an options/parametrization structure
            EXT2 = EXT2(logical(ModelParameter{kk,1}.BasisSets.includeInFit(ss,:)));
            parametrization.lorentzLB.ex    = EXT2;
            parametrization.lorentzLB.init  = EXT2;
            SDT2 = SDT2(logical(ModelParameter{kk,1}.BasisSets.includeInFit(ss,:)));
            parametrization.lorentzLB.sd = SDT2;
            SDSH = SDSH(logical(ModelParameter{kk,1}.BasisSets.includeInFit(ss,:)));
            parametrization.freqShift.sd = SDSH;
        else
            % ModelParameter{kk,1}.indexMMLipBasisFunctionInBasis;
            % ModelParameter{kk,1}.BasisSets.indexMMLipBasisFunction = zeros(size(logical(ModelParameter{kk,1}.BasisSets.includeInFit(ss,:))));
            % Only initialize EX/SD values (for frequency shift and
            % Lorentzian LB) for the metabolites:
            [EXT2, SDT2, SDSH] = fit_setExSDValues(DataToModel{kk}, basisSet);
             % Store EXT2, SDT2, SDSH to an options/parametrization structure
            EXT2 = EXT2(logical(ModelParameter{kk,1}.BasisSets.includeInFit(ss,:)) & logical(~ModelParameter{kk,1}.BasisSets.indexMMLipBasisFunction));
            parametrization.lorentzLB.ex    = EXT2;
            parametrization.lorentzLB.init  = EXT2;
            SDT2 = SDT2(logical(ModelParameter{kk,1}.BasisSets.includeInFit(ss,:)) & logical(~ModelParameter{kk,1}.BasisSets.indexMMLipBasisFunction));
            parametrization.lorentzLB.sd = SDT2;
            SDSH = SDSH(logical(ModelParameter{kk,1}.BasisSets.includeInFit(ss,:)) & logical(~ModelParameter{kk,1}.BasisSets.indexMMLipBasisFunction));
            parametrization.freqShift.sd = SDSH;
        end
       
        % Update the FitObject
        ModelParameter{kk,1}.updateParametrization(parametrization);

        % Is this modeling or optimization of the regularization parameter
        if ~strcmp(ModelProcedure.Steps{ss}.module,'OptimReg')
            fprintf('Running model of dataset #%i of %i \n', kk, length(DataToModel));
            % Run steps defined ModelProcedure struct
            % Loop across all steps with anonymous calls defined in the ModelProcedure struct
            ModelParameter{kk}.createModel;
        else
            fprintf('Optimize regularizer of dataset #%i of %i \n', kk, length(DataToModel));
            % Run steps defined in ModelProcedure struct N times to
            % optimize the regularization parameter
            opts.regularizer    = ModelProcedure.Steps{ss}.regularizer;
            if ss == 1
                ModelParameter{kk,1}.optimizeRegularization(opts);
            else
                ModelParameter{kk,1}.updateOptsAccordingToStep(opts);
                ModelParameter{kk,1}.optimizeRegularization(opts);
            end
        end
    end
    if DoDataScaling >= 1                                                   % We did data scaling - so let's store the scale
        ModelParameter{kk,1}.scale = scaleData(kk);
    end
    DataToModel{kk} 	= [];                                               % Memory efficiency
end

%% 4. Save parameter results (and export as NII)
% Save ModelParameter to be includeded in the container and final
% results in NII LCM format for easy reading and visualization

% Lets reduce the file size
for kk = 2 : size(ModelParameter,1)
    ModelParameter{kk,1}.economizeStorage(1,1);                         % Remove basis set and jacobians
end
ModelParameter{1,1}.economizeStorage(0,1);                              % Keep first basis set but no jacobians

end

%% Helper functions below

function isStructOrCell(f)
assert(iscell(f) || isstruct(f));
end

function isCharOrStruct(f)
assert(ischar(f) || isstruct(f));
end

function fwhm = measureFWHM(ppm, fid)
% Measures FWHM of singlet signal
spec = fftshift(fft(fid));
[maxAmp, maxIdx] = max(real(spec));
halfMax = maxAmp/2;
specShifted = abs(real(spec)-halfMax);
[~,hmIdx] = min(specShifted);

fwhm = 2*abs(ppm(maxIdx) - ppm(hmIdx));
end
