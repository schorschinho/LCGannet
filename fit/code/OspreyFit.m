function [MRSCont] = OspreyFit(MRSCont)
%% [MRSCont] = OspreyFit(MRSCont)
%   This function performs spectral fitting on MRS data loaded previously
%   using OspreyLoad.
%
%   The method of fit, fitting range, included metabolites and other
%   settings are set in the job file.
%
%   USAGE:
%       MRSCont = OspreyFit(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-24)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-24: First version of the code.

outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));

if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

% ----- Load fit settings and fit the metabolite data -----
% Checking for version, toolbox, and previously run modules
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyFit');
% Start timer
MRSCont.runtime.Fit = 0;

% Initialise the fit - this step includes:
% - Parse the correct basis set
% - Apply settings on which metabolites/MM/lipids to include in the fit
% - Check for inconsistencies between basis set and data
[MRSCont] = osp_fitInitialise(MRSCont);


% New Osprey gLCM model
if strcmpi(MRSCont.opts.fit.method, 'Osprey_gLCM')
    [MRSCont] = osp_fit_gLCM(MRSCont);
else
    % Call the fit functions (depending on sequence type)
    if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
        if MRSCont.flags.isUnEdited
            [MRSCont] = osp_fitUnEdited(MRSCont);
        elseif MRSCont.flags.isMEGA
            [MRSCont] = osp_fitMEGA(MRSCont);
        elseif MRSCont.flags.isHERMES
            [MRSCont] = osp_fitHERMES(MRSCont);
        elseif MRSCont.flags.isHERCULES
            % For now, fit HERCULES like HERMES data
            [MRSCont] = osp_fitHERCULES(MRSCont);
        else
            msg = 'No flag set for sequence type!';
            fprintf(msg);
            error(msg);
        end
    else
        [MRSCont] = osp_fitMultiVoxel(MRSCont);
    end
end

% ----- Perform water reference and short-TE water fit -----
% The water signal is automatically integrated when the LCModel fit option is
% being used. In Osprey, we explicitly model the water data with a
% dedicated simulated water basis function.
if strcmpi(MRSCont.opts.fit.method, 'Osprey')

    % If water reference exists, fit it
    if MRSCont.flags.hasRef
        refFitTime = tic;
        % Loop over all the datasets here
        for kk = 1:MRSCont.nDatasets(1)
            [~] = printLog('OspreyFitRef',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
            if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.ref.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
                [MRSCont] = osp_fitWater(MRSCont, kk, 'ref');
            end
        end
        time = toc(refFitTime);
        if MRSCont.flags.isGUI
            set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
            pause(1);
        end
        fprintf('... done.\n Elapsed time %f seconds\n',time);
        MRSCont.runtime.FitRef = time;
        MRSCont.runtime.Fit = MRSCont.runtime.Fit + time;
    end

    % If short TE water reference exists, fit it
    if MRSCont.flags.hasWater
        waterFitTime = tic;
        % Loop over all the datasets here
        for kk = 1:MRSCont.nDatasets(1)
            [~] = printLog('OspreyFitWater',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
            if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.w.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
                [MRSCont] = osp_fitWater(MRSCont, kk, 'w');
            end
        end
        time = toc(waterFitTime);
        fprintf('... done.\n Elapsed time %f seconds\n',time);
        MRSCont.runtime.FitWater = time;
        MRSCont.runtime.Fit = MRSCont.runtime.Fit + time;
    end

    MRSCont.runtime.Fit = MRSCont.runtime.Fit + MRSCont.runtime.FitMet;

else if strcmpi(MRSCont.opts.fit.method, 'Osprey_gLCM')
  
        % If water reference exists, fit it
        if MRSCont.flags.hasRef
           % We want a loop over the extra dimension for separate fitting
            SeparateExtraDims = 1;
            if MRSCont.processed.metab{1}.dims.extras > 0
                SeparateExtraDims = MRSCont.processed.metab{1}.sz(MRSCont.processed.ref{1}.dims.extras);
            end
            % Read model procedure 
            ModelProcedure = jsonToStruct(MRSCont.opts.fit.ModelProcedure.ref{1});
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end
            if ModelProcedure.Steps{1, 1}.extra.flag                                    % Model is actually 2D along extra dimension so we don't want to loop 
                SeparateExtraDims = 1;
            end 
            refFitTime = tic;
            for ex = 1 : SeparateExtraDims          
                if ~isfield(ModelProcedure,'basisset') || ~isfield(ModelProcedure.basisset, 'file') || ... 
                    isempty(ModelProcedure.basisset.file)
                    if ~iscell(MRSCont.fit.basisSet)
                        ModelProcedure.basisset.file = {MRSCont.fit.basisSet};
                    else
                        ModelProcedure.basisset.file = MRSCont.fit.basisSet;
                    end
                end
                if SeparateExtraDims > 1
                    ModelProcedure.basisset.opts.index = ex;
                end
                % Loop over all the datasets here
                if isprop(MRSCont.fit.results.metab{1,1}, 'scale')
                    scale = [];
                    for kk = 1:MRSCont.nDatasets(1)
                        scale = [scale MRSCont.fit.results.metab{1,kk,1,ex}.scale];
                    end
                else
                    scale = 0;
                end      
                [MRSCont.fit.results.ref(1,:,1,ex)] = Osprey_gLCM(MRSCont.processed.ref,ModelProcedure,0,0,scale)';         
            end
        time = toc(refFitTime);
        if MRSCont.flags.isGUI
            set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
            pause(1);
        end
        fprintf('... done.\n Elapsed time %f seconds\n',time);
        MRSCont.runtime.FitRef = time;
        MRSCont.runtime.Fit = MRSCont.runtime.Fit + time;
    end

    % If short TE water reference exists, fit it
    if MRSCont.flags.hasWater
        % We want a loop over the extra dimension for separate fitting
        SeparateExtraDims = 1;
        if MRSCont.processed.metab{1}.dims.extras > 0
            SeparateExtraDims = MRSCont.processed.w{1}.sz(MRSCont.processed.metab{1}.dims.extras);
        end
        % Read model procedure 
        ModelProcedure = jsonToStruct(MRSCont.opts.fit.ModelProcedure.w{1});
        if isstruct(ModelProcedure.Steps)
            ModelProcedureCell = cell(size(ModelProcedure.Steps));
            for ss = 1 : size(ModelProcedure.Steps,1)
                ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
            end
            ModelProcedure.Steps = ModelProcedureCell;
        end
        if ModelProcedure.Steps{1, 1}.extra.flag                                    % Model is actually 2D along extra dimension so we don't want to loop 
            SeparateExtraDims = 1;
        end 
        waterFitTime = tic;
        for ex = 1 : SeparateExtraDims
            if ~isfield(ModelProcedure,'basisset') || ~isfield(ModelProcedure.basisset, 'file') || ... 
                isempty(ModelProcedure.basisset.file)
                if ~iscell(MRSCont.fit.basisSet)
                    ModelProcedure.basisset.file = {MRSCont.fit.basisSet};
                else
                    ModelProcedure.basisset.file = MRSCont.fit.basisSet;
                end
            end
            if SeparateExtraDims > 1
                ModelProcedure.basisset.opts.index = ex;
            end
            % Loop over all the datasets here
            if isprop(MRSCont.fit.results.metab{1,1}, 'scale')
                scale = [];
                for kk = 1:MRSCont.nDatasets(1)
                    scale = [scale MRSCont.fit.results.metab{1,kk,1,ex}.scale];
                end
            else
                scale = 0;
            end      
            [MRSCont.fit.results.w(1,:,1,ex)] = Osprey_gLCM(MRSCont.processed.w,ModelProcedure,0,0,scale)';         
        end
        time = toc(waterFitTime);
        fprintf('... done.\n Elapsed time %f seconds\n',time);
        MRSCont.runtime.FitWater = time;
        MRSCont.runtime.Fit = MRSCont.runtime.Fit + time;
    end

    MRSCont.runtime.Fit = MRSCont.runtime.Fit + MRSCont.runtime.FitMet;

    else
        MRSCont.runtime.Fit =  MRSCont.runtime.FitMet;
    end
end
[~] = printLog('Fulldone',MRSCont.runtime.Fit,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);


%% If DualVoxel or MRSI we want to extract y-axis scaling
if MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI
    MRSCont = osp_scale_yaxis(MRSCont,'OspreyLoad');
    MRSCont.fit.resBasisSet = MRSCont.fit.resBasisSet{2,2};
end

%% Delete redundant resBasiset entries
if strcmpi(MRSCont.opts.fit.method, 'Osprey')
    if ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        FitNames = fieldnames(MRSCont.fit.resBasisSet);
        NoFit = length(fieldnames(MRSCont.fit.resBasisSet));
        for sf = 1 : NoFit
            MRSCont.fit.resBasisSet.(FitNames{sf}) = MRSCont.fit.resBasisSet.(FitNames{sf})(:,MRSCont.info.(FitNames{sf}).unique_ndatapoint_spectralwidth_ind,:);
            for combs = 1 : length(MRSCont.info.(FitNames{sf}).unique_ndatapoint_spectralwidth_ind)
                resBasisSetNew.(FitNames{sf}).([MRSCont.info.(FitNames{sf}).unique_ndatapoint_spectralwidth{combs}]) = MRSCont.fit.resBasisSet.(FitNames{sf})(:,combs,:);
            end
        end
        MRSCont.fit.resBasisSet = resBasisSetNew;
    end
end

%% Store  and print some QM parameters
if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_fit_Quality(MRSCont);

    L = length(MRSCont.QM.tables.Properties.VariableNames);
    % Store data quality measures in csv file
    if MRSCont.flags.isUnEdited
        relResA = MRSCont.QM.relAmpl.metab_A';
        MRSCont.QM.tables.relResA = relResA;
    elseif MRSCont.flags.isMEGA
        if strcmp( MRSCont.opts.fit.style, 'Separate')
            relResA = MRSCont.QM.relAmpl.metab_A';
            relResdiff1 = MRSCont.QM.relAmpl.metab_diff1';
            MRSCont.QM.tables.relResA = relResA;
            MRSCont.QM.tables.relResdiff1 = relResdiff1;
        else
            relRessum = MRSCont.QM.relAmpl.metab_sum';
            relResdiff1 = MRSCont.QM.relAmpl.metab_diff1';
            MRSCont.QM.tables.relRessum = relRessum;
            MRSCont.QM.tables.relResdiff1 = relResdiff1;
        end
    elseif MRSCont.flags.isHERMES
            relRessum = MRSCont.QM.relAmpl.metab_sum';
            relResdiff1 = MRSCont.QM.relAmpl.metab_diff1';
            relResdiff2 = MRSCont.QM.relAmpl.metab_diff2';
            MRSCont.QM.tables.relRessum = relRessum;
            MRSCont.QM.tables.relResdiff1 = relResdiff1;
            MRSCont.QM.tables.relResdiff2 = relResdiff2;
    elseif MRSCont.flags.isHERCULES
        % For now, process HERCULES like HERMES data
            relRessum = MRSCont.QM.relAmpl.metab_sum';
            relResdiff1 = MRSCont.QM.relAmpl.metab_diff1';
            relResdiff2 = MRSCont.QM.relAmpl.metab_diff2';
            MRSCont.QM.tables.relRessum = relRessum;
            MRSCont.QM.tables.relResdiff1 = relResdiff1;
            MRSCont.QM.tables.relResdiff2 = relResdiff2;
    else
        msg = 'No flag set for sequence type!';
        fprintf(msg);
        error(msg);
    end

    % Loop over field names to populate descriptive fields of table for JSON export
    for JJ = L:length(MRSCont.QM.tables.Properties.VariableNames)
        switch MRSCont.QM.tables.Properties.VariableNames{JJ}
            case 'relResA'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'relResA'} = 'relResA';%CWDJ??
                MRSCont.QM.tables.Properties.VariableDescriptions{'relResA'} = 'Fit quality number for spectrum A relative amplitude of the residual compared to the standard deveiation of the noise';
                MRSCont.QM.tables.Properties.VariableUnits{'relResA'} = 'arbitrary';
            case 'relRessum'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'relRessum'} = 'relRessum';
                MRSCont.QM.tables.Properties.VariableDescriptions{'relRessum'} = 'Fit quality number for sum spectrum relative amplitude of the residual compared to the standard deveiation of the noise';
                MRSCont.QM.tables.Properties.VariableUnits{'relRessum'} = 'arbitrary';
            case 'relResdiff1'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'relResdiff1'} = 'relResdiff1';
                MRSCont.QM.tables.Properties.VariableDescriptions{'relResdiff1'} = 'Fit quality number for diff1 spectrum relative amplitude of the residual compared to the standard deveiation of the noise';
                MRSCont.QM.tables.Properties.VariableUnits{'relResdiff1'} = 'arbitrary';
            case 'relResdiff2'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'relResdiff2'} = 'relResdiff2';
                MRSCont.QM.tables.Properties.VariableDescriptions{'relResdiff2'} = 'Fit quality number for diff2 spectrum relative amplitude of the residual compared to the standard deveiation of the noise';
                MRSCont.QM.tables.Properties.VariableUnits{'relResdiff2'} = 'arbitrary';
        end
    end

    %Output as .tsv
    osp_WriteBIDsTable(MRSCont.QM.tables, [outputFolder filesep 'QM_processed_spectra'])
end

%% Clean up and save
% Set exit flags and version
MRSCont.flags.didFit           = 1;

diary off
% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyFit');
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end
