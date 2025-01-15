function out = osp_plotFit(MRSCont, kk, which_spec,VoxelIndex, conc, stagFlag, xlab, ylab, figTitle)
%% out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%   Creates a figure showing data stored in an Osprey data container, as
%   well as the fit to it, the baseline, the residual, and contributions
%   from the individual metabolites.
%
%   USAGE:
%       out = osp_plotFit(MRSCont, kk, which, GUI, conc, stagFlag, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'metab' (default)
%                               'ref'
%                               'w'
%                               'mm' re_mm
%       VoxelIndex = Vector with three entries experiment number (1), subspectrum (2),
%                    and basisset (3).
%                   OPTIONS: experiment number only supports 1 for now
%                            subspectra = 1 for unedited, off (MEGA), sum (HERMES/HERCULES)
%                            subspectra = 2 for diff1 (MEGA), diff2 (HERMES/HERCULES)
%                            subspectra = 3 for diff3 (HERMES/HERCULES)
%                            basisset = 1 for default parameterization
%                            basisset = 2 for subject-specific measured MM
%                            (if metabolite-nulled spectra are supplied)
%                            basisset = 3 for cohort-mean measured MM
%                            (if metabolite-nulled spectra are supplied & opts.fit.MeanMM = 1)
%                            
%       conc      = flag if concatenate fitting was used
%       stagFlag  = flag to decide whether basis functions should be plotted
%                   vertically staggered or simply over one another
%                   (optional)
%                    OPTIONS:   1 = staggered (default)
%                               0 = not staggered
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    error('Trying to plot fitted data, but fit has not been performed yet. Run OspreyFit first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided
if nargin<9
    [~,filen,ext] = fileparts(MRSCont.files{kk});
    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        if strcmp(which_spec, 'conc')
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' conc ' fit plot:\n' filen ext]);
        else
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' which_spec ' fit plot:\n' filen ext]);
        end
    elseif  (MRSCont.flags.isPRIAM && isfield(MRSCont.flags,'isPRIAM'))
        if nargin<4
            VoxelIndex = 1;
        end
        if nargin<5
            conc = 'diff1';
        end
        if strcmp(which_spec, 'conc')
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' conc ' fit plot:\n' filen ext '\n Voxel ' num2str(VoxelIndex)]);
        else
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' which_spec ' fit plot:\n' filen ext  '\n Voxel ' num2str(VoxelIndex)]);
        end
    else
        if nargin<4
            VoxelIndex = [1 1 1];
        end
        if strcmp(which_spec, 'conc')
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' conc ' fit plot:\n' filen ext '\n Voxel ' num2str(VoxelIndex(1)) ' ' num2str(VoxelIndex(2))]);
        else
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' which_spec ' fit plot:\n' filen ext  '\n Voxel ' num2str(VoxelIndex(1)) ' ' num2str(VoxelIndex(2))]);
        end
    end
    if nargin<8
        ylab='';
        if nargin<7
            xlab='chemical shift (ppm)';
            if nargin<6
                stagFlag = 1;
                if nargin<5
                    conc = 'diff1';
                    if nargin<4
                        VoxelIndex = [1 1 1];
                        if nargin < 3
                            which_spec = 'off';
                            if nargin < 2
                                kk = 1;
                                if nargin<1
                                    error('ERROR: no input Osprey container specified.  Aborting!!');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


%%% 2. EXTRACT DATA TO PLOT %%%
% We will do this very differently for Osprey and LCModel:
% For Osprey, we will take the basis set and apply the model parameters to
% produce plots on the fly.
% For LCModel, we will load the model output that is provided to us through
% the .coord and .table files.

switch fitMethod
    case 'Osprey'
        % Extract processed spectra and fit parameters
        if (MRSCont.flags.isPRIAM == 1) || (MRSCont.flags.isMRSI == 1)
            if ~exist('VoxelIndex') && (MRSCont.flags.isPRIAM == 1)
                VoxelIndex = 1;
            elseif ~exist('VoxelIndex') && (MRSCont.flags.isMRSI == 1)
                VoxelIndex = [1 1 1];
            end
            if  strcmp(which_spec, 'conc')
                dataToPlot=op_takeVoxel(MRSCont.processed.(conc){kk},VoxelIndex);
            else
                if strcmp(which_spec, 'off')
                    dataToPlot=op_takeVoxel(MRSCont.processed.A{kk},VoxelIndex);
                else
                    dataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);
                end
            end

            if (MRSCont.flags.isPRIAM == 1)
                if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
                    fitRangePPM = MRSCont.opts.fit.rangeWater;
                    basisSet    = MRSCont.fit.resBasisSet{VoxelIndex}.(which_spec).water.(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]);
                else if strcmp(which_spec, 'conc')
                        fitRangePPM = MRSCont.opts.fit.range;
                        basisSet    = MRSCont.fit.resBasisSet{VoxelIndex}.(which_spec).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]);
                    else
                        fitRangePPM = MRSCont.opts.fit.range;
                        basisSet    = MRSCont.fit.resBasisSet{VoxelIndex}.(which_spec).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]);
                    end
                end
            else
                if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
                    fitRangePPM = MRSCont.opts.fit.rangeWater;
                    basisSet    = MRSCont.fit.resBasisSet.(which_spec).water.(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]);
                else if strcmp(which_spec, 'conc')
                        fitRangePPM = MRSCont.opts.fit.range;
                        basisSet    = MRSCont.fit.resBasisSet.(which_spec).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]);
                    else
                        fitRangePPM = MRSCont.opts.fit.range;
                        basisSet    = MRSCont.fit.resBasisSet.(which_spec).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]);
                    end
                end
            end
        else
            if  strcmp(which_spec, 'conc')
                dataToPlot  = MRSCont.processed.(conc){kk};
            else
                dataToPlot = MRSCont.processed.(which_spec){kk};
            end

            if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
                fitRangePPM = MRSCont.opts.fit.rangeWater;
                basisSet    = MRSCont.fit.resBasisSet.(which_spec).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){VoxelIndex,1};
            else if strcmp(which_spec, 'conc')
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet.(which_spec).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){VoxelIndex,1};
                else
                    fitRangePPM = MRSCont.opts.fit.range;
                    
                    if strcmp(which_spec,'mm') && MRSCont.flags.isUnEdited
                        fitRangePPM = [0.2 4.2];
                    end
                    basisSet    = MRSCont.fit.resBasisSet.([which_spec]).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){VoxelIndex(3),1,VoxelIndex(2)};
                    if VoxelIndex(3) == 2 % We need to insert the subject specific MM basis function into the basis set
                        if VoxelIndex(2)==1
                            index = find(strcmp(MRSCont.processed.mm{kk}.names,'A_spline')); 
                        end
                        if VoxelIndex(2)==2
                            index = find(strcmp(MRSCont.processed.mm{kk}.names,'diff1_spline')); 
                        end
                        mm_clean_spline = op_takesubspec(MRSCont.processed.mm{kk},index);
                        mm_clean_spline               = op_zeropad(mm_clean_spline, 2);  
                        ind = find(strcmp(basisSet.name,'MMExp'));
                        basisSetfactor = op_freqrange(basisSet,0,1.2);
                        mm_clean_spline_factor = op_freqrange(mm_clean_spline,0.7,1.1);
                        factor = (max(real(basisSetfactor.specs(:,ind)))/max(real(mm_clean_spline_factor.specs)));
                        mm_clean_spline = op_ampScale(mm_clean_spline,factor);
                        basisSet.fids(:,ind) = mm_clean_spline.fids;
                        basisSet.specs(:,ind) = mm_clean_spline.specs;
                    end
                end
            end
            if  strcmp(which_spec, 'conc')
                dataToPlot  = MRSCont.processed.(conc){kk};
            else
                dataToPlot   = op_takesubspec(MRSCont.processed.(which_spec){kk},find(strcmp(MRSCont.processed.(which_spec){kk}.names,basisSet.names{1})));
            end
        end



        % Get the fit parameters

        if (MRSCont.flags.isPRIAM == 1)
            fitParams   = MRSCont.fit.results{VoxelIndex}.(which_spec).fitParams{kk};
        elseif (MRSCont.flags.isMRSI == 1)
            fitParams   = MRSCont.fit.results{VoxelIndex(1), VoxelIndex(2)}.(which_spec).fitParams{kk};
        elseif strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            fitParams   = MRSCont.fit.results.(which_spec).fitParams{VoxelIndex(3),kk};
        else
            fitParams   = MRSCont.fit.results.(which_spec).fitParams{VoxelIndex(3),kk,VoxelIndex(2)};
        end
        % Pack up into structs to feed into the reconstruction functions
        inputData.dataToFit                 = dataToPlot;
        inputData.basisSet                  = basisSet;
%         if (length(fitParams.ampl) == 3)
%             inputData.basisSet_mm                  = MRSCont.fit.basisSet_mm;
%         end
        if (MRSCont.flags.isPRIAM == 1)
            inputSettings.scale                 = MRSCont.fit.scale{kk};
        else
            inputSettings.scale                 = MRSCont.fit.scale{kk};
        end
        inputSettings.fitRangePPM           = fitRangePPM;
        inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
        inputSettings.fitStyle              = MRSCont.opts.fit.style;
        inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
        inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
        inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
        inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;
        inputSettings.concatenated.Subspec  = conc;
        if isfield(MRSCont.opts.fit,'GAP') && ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
            inputSettings.GAP = MRSCont.opts.fit.GAP.(dataToPlot.names{1});
        else
            inputSettings.GAP = [];
        end

    case 'LCModel'
        fitRangePPM = MRSCont.opts.fit.range;

        % We do not store separate fit parameters for the water fit;
        % instead, we have saved the unsuppressed water area in the 'off'
        % parameter set.
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')

            % Do nothing for now. We'll load the water spectrum in the next
            % step.

        else
            % Get the fit parameters; they will be passed to the model
            % output in the next step
            if (MRSCont.flags.isPRIAM == 1)
                fitParams   = MRSCont.fit.results{VoxelIndex}.(which_spec).fitParams{kk};
            elseif (MRSCont.flags.isMRSI == 1)
                fitParams   = MRSCont.fit.results{VoxelIndex(1), VoxelIndex(2)}.(which_spec).fitParams{kk};
            else
                fitParams   = MRSCont.fit.results.(which_spec).fitParams{VoxelIndex(3),kk,VoxelIndex(2)};
            end

        end

end

if isfield(MRSCont.opts, 'cosmetics')
    fitParams.lb = MRSCont.opts.cosmetics.LB;
else
    fitParams.lb = 0;
    MRSCont.opts.cosmetics.Zoom = 1;
end

%%% 3. PREPARE LINES TO DISPLAY %%%
% Extract data, ppm axes, fit, residual, baseline, and individual
% metabolite contributions.
% Obviously, depending on the fit method that has been used, we need to
% reconstruct the plots with fit method functions
switch fitMethod
    % Depending on whether metabolite or water data are to be
    % displayed, create the plots via different models
    case 'Osprey'
        
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            if strcmp(inputSettings.fitStyle,'Concatenated')
                [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
            else
                [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
            end
        end
    case 'OspreyAsym'
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            if strcmp(inputSettings.fitStyle,'Concatenated')
                [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
            else
                [ModelOutput] = fit_OspreyAsymParamsToModel(inputData, inputSettings, fitParams);
            end
        end
    case 'OspreyNoLS'
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            if strcmp(inputSettings.fitStyle,'Concatenated')
                [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
            else
                [ModelOutput] = fit_OspreyNoLSParamsToModel(inputData, inputSettings, fitParams);
            end
        end
    case 'LCModel'
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            % Just load the water spectrum (we don't have a fit)
            waterSpec         = op_freqrange(MRSCont.processed.(which_spec){kk}, 0, 9);
            ModelOutput.data  = real(waterSpec.specs);
            ModelOutput.ppm   = waterSpec.ppm;
        else
            % Get the LCModel plots we previously extracted from .coord
            % etc.
            [ModelOutput] = fit_LCModelParamsToModel(fitParams);
        end

end

%re_mm
%For MM, prepare a 'clean' MM spectrum that has metabolite signals pulled
%out.
if (strcmp(which_spec, 'mm'))
    Met_corr_spectrum  = sum(ModelOutput.indivMets(:,1:basisSet.nMets),2);
end
if (strcmp(which_spec, 'diff1_mm'))
   Met_corr_spectrum  = sum(ModelOutput.indivMets(:,1:basisSet.nMets),2);
end


%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
canvasSize  = get(0,'defaultfigureposition');
if stagFlag && ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
    canvasSize(4) = canvasSize(4) * 1.8;
end
out = figure('Position', canvasSize);
% Prepare a couple of useful variables
switch fitMethod
    case 'Osprey'
        nBasisFct = basisSet.nMets;
        if isfield(basisSet, 'nMM')
            nBasisFct = nBasisFct + basisSet.nMM;
        end
    case 'LCModel'
        % Number of metabolites and lipid/MM basis functions
        nBasisFct = length(fitParams.name);
        nComb   = sum(~cellfun(@isempty, strfind(fitParams.name, '_')));
        nBasisFct =  nBasisFct - nComb; % We don't plot the combinations
end


%%% 5. PLOT DATA, FIT, RESIDUAL, BASELINE %%%
% Create a fitMethod-specific theme
switch fitMethod
    case 'Osprey'
        colorData = MRSCont.colormap.Foreground;
        colorFit  = MRSCont.colormap.Accent;
        linewidthFit = 1.6;
        linewidthResidual = 1;
        colorBaseline = MRSCont.colormap.LightAccent;
    case 'LCModel'
        colorData = MRSCont.colormap.Foreground;
        colorFit  = 'r';
        linewidthFit = 1.2;
        linewidthResidual = 0.5;
        colorBaseline = MRSCont.colormap.Foreground;
end

% Unpack the ModelOutput struct
ppm         = ModelOutput.ppm;
dataToPlot  = ModelOutput.data;

switch fitMethod
    case 'Osprey'
        % If water, don't get baseline and individual fits
        fit         = ModelOutput.completeFit;
        residual    = ModelOutput.residual;
        if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
            baseline    = ModelOutput.baseline;
            indivPlots  = ModelOutput.indivMets;
            basisSetNames = basisSet.name;
        end
    case 'LCModel'
        % If water, don't get anything except the data
        if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
            fit         = ModelOutput.completeFit;
            residual    = ModelOutput.residual;
            baseline    = ModelOutput.baseline;
            indivPlots  = ModelOutput.indivMets;
            basisSetNames = fitParams.name;
            indivPlots(:,contains(basisSetNames,'_')) = [];
            basisSetNames(contains(basisSetNames,'_')) = [];
        end

end


if isfield(MRSCont.plot,'fit') && MRSCont.plot.fit.match
    if strcmp(which_spec, 'conc')
        stagData = MRSCont.plot.fit.(conc).stagData(kk);
        maxPlot = MRSCont.plot.fit.(conc).maxPlot(kk);
    else
        stagData = MRSCont.plot.fit.(which_spec).stagData(kk);
        maxPlot = MRSCont.plot.fit.(which_spec).maxPlot(kk);
    end
else
    % Determine a positive stagger to offset data, fit, residual, and
    % baseline from the individual metabolite contributions
    stagData = 0.1*(max(abs(min(dataToPlot)), abs(max(dataToPlot))))/MRSCont.opts.cosmetics.Zoom;
    maxPlot = max(dataToPlot + abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData;
end
% Add the data and plot
hold on;
plot(ppm, (zeros(1,length(ppm)) + stagData)/maxPlot, 'Color', colorData); % Zeroline
if (contains(which_spec, 'mm'))
   plot(ppm, (dataToPlot + stagData-Met_corr_spectrum)/maxPlot, 'Color',[1 0 0.1]); % Data
end
plot(ppm, (dataToPlot + stagData)/maxPlot, 'Color', colorData); % Data
plot(ppm, (fit + stagData)/maxPlot, 'Color', colorFit, 'LineWidth', linewidthFit); % Fit
plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot) + stagData)/maxPlot, 'Color', colorData, 'LineWidth', 1); % Maximum Data

plot(ppm, (residual + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color', colorData, 'LineWidth', linewidthResidual); % Residual
plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color', colorData, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot, 'Color', colorData, 'LineWidth', 1); % Max Residue



text(fitRangePPM(1), (0 + stagData)/maxPlot, '0', 'FontSize', 10,'Color', MRSCont.colormap.Foreground); %Zeroline text
text(fitRangePPM(1), (0 + max(dataToPlot) + stagData)/maxPlot-0.05, num2str(max(dataToPlot),'%1.2e'), 'FontSize', 10,'Color', MRSCont.colormap.Foreground); % Maximum Data Text
text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, '0', 'FontSize', 10,'Color', MRSCont.colormap.Foreground); %Zeroline Residual text
text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot +0.05, num2str(abs(max(dataToPlot - fit)),'%1.2e'), 'FontSize', 10,'Color', MRSCont.colormap.Foreground); %Max Residue text

if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w') || contains(which_spec, 'mm')) 
    plot(ppm, (real(baseline) + stagData)/maxPlot, 'k', 'LineWidth', 1, 'Color', colorBaseline);
end


%%% 6. PLOT BASIS FUNCTIONS %%%
% Plot separate metabolite basis functions only if not water
if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
    if stagFlag
        % Staggered plots will be in all black and separated by the mean of the
        % maximum across all spectra
        %         stag = max(abs(mean(max(real(appliedBasisSet.specs)))), abs(mean(min(real(appliedBasisSet.specs))))) * MRSCont.fit.scale{kk};
        stag = maxPlot *  (3.5 - MRSCont.opts.cosmetics.Zoom) / nBasisFct;
        % Loop over all basis functions

        for rr = 1:nBasisFct
            % Instead of a MATLAB legend, annotate each line separately with the
            % name of the metabolite
            plot(ppm, (indivPlots(:,rr) - rr*stag)/maxPlot, 'Color', MRSCont.colormap.Foreground);
            text(fitRangePPM(1), (- rr*stag)/maxPlot, basisSetNames{rr}, 'FontSize', 10,'Color', MRSCont.colormap.Foreground);
        end

        % Preliminary formatting; might need some more stability here, or
        % differentiation based on sequence type
        set(gca, 'YLim', [(-nBasisFct-1)*stag/maxPlot  1]);
        hold off;

    else
        % If not staggered, plots will simply be made with distinguishable
        % colors
        colours = distinguishable_colors(nBasisFct);

        % Loop over all basis functions
        for rr = 1:nBasisFct
            plot(ppm, indivPlots(:,rr)/maxPlot, 'Color', colours(rr,:), 'LineWidth', 1);
        end
        legend([newline, basisSetNames], 'Orientation', 'horizontal');
        hold off

    end

else
    % Preliminary formatting; might need some more stability here, or
    % differentiation based on sequence type
    set(gca, 'YLim', [0  1.2]);
    hold off;
    % If water is being shown, show a simple legend
    %legend('Data', 'Fit', 'Residual');

end


%%% 7. DESIGN FINETUNING %%%
% Adapt common style for all axes
set(gca, 'XDir', 'reverse', 'XLim', [fitRangePPM(1), fitRangePPM(end)], 'XMinorTick', 'On');
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
ticks = get(gca,'XTick');
set(gca, 'XTick', unique(round(ticks)));
% If no y caption, remove y axis
if isempty(ylab)
    if ~MRSCont.flags.isGUI
        set(gca, 'YColor', 'w');
        % Black axes, white background
        set(gca, 'XColor', 'k');
        set(gca, 'Color', 'w');
        set(gcf, 'Color', 'w');
        title(figTitle, 'Interpreter', 'none');
    else
        set(gca, 'YColor', MRSCont.colormap.Background);
        set(gca,'YTickLabel',{});
        set(gca,'YTick',{});
        % Dirtywhite axes, light gray background
        set(gca, 'XColor', MRSCont.colormap.Foreground);
        set(gca, 'Color', MRSCont.colormap.Background);
        set(gcf, 'Color', MRSCont.colormap.Background);
        title(figTitle, 'Interpreter', 'none', 'Color', MRSCont.colormap.Foreground);
    end
else
    set(gca, 'YColor', 'k');
end

box off;
xlabel(xlab, 'FontSize', 16);
ylabel(ylab, 'FontSize', 16);


%%% 8. ADD OSPREY LOGO AND TIGHTEN FIGURE %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
else
    out = tightfig(out);
end

end