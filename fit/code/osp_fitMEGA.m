function [MRSCont] = osp_fitMEGA(MRSCont)
%% [MRSCont] = osp_fitMEGA(MRSCont)
%   This function performs spectral fitting of MEGA-edited MRS data.
%
%   USAGE:
%       [MRSCont] = osp_fitMEGA(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-12)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-04-12: First version of the code.


% Loop over all the datasets here
metFitTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

if MRSCont.flags.hasMM == 1
    % Load basis set, calculate specs, and sort:
    basisSet_mm = load(MRSCont.opts.fit.basisSetFile);
    basisSet_mm = basisSet_mm.BASIS;
    basisSet_mm = osp_recalculate_basis_specs(basisSet_mm);
    basisSet_mm = fit_sortBasisSet(basisSet_mm);
    
    % Extract edit-off spectrum from the basis:
    metabList_mm = fit_createMetabListMM('unedited');
    basisSet_mm_A = fit_selectMetabs(basisSet_mm, metabList_mm, 1);
    basisSet_mm_A.fids = basisSet_mm_A.fids(:,:,1);
    basisSet_mm_A.specs = basisSet_mm_A.specs(:,:,1);
    
    % Extract difference spectrum from the basis:
    metabList_mm = fit_createMetabListMM('MEGA');
    basisSet_mm_diff1 = fit_selectMetabs(basisSet_mm, metabList_mm, 1);
    basisSet_mm_diff1.fids = basisSet_mm_diff1.fids(:,:,3);
    basisSet_mm_diff1.specs = basisSet_mm_diff1.specs(:,:,3);
end

for kk = 1:MRSCont.nDatasets(1)
    % ----- Osprey fit pipeline -----
    if strcmpi(MRSCont.opts.fit.method, 'Osprey')
        [~] = printLog('OspreyFit',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
    
        %%% 1. DETERMINE THE FITTING STYLE %%%
        % Extract fit options
        fitOpts     = MRSCont.opts.fit;
        fitModel    = fitOpts.method;
        fitStyle    = fitOpts.style;
    
        %%% 2. SEPARATE FIT %%%
        % For the separate (classic) MEGA fit, model the EDIT-OFF and DIFF
        % spectra separately.
        if strcmp(fitStyle, 'Separate')
            if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > size(MRSCont.fit.results.metab.fitParams,2))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
                %%% 2a. FIT OFF-SPECTRUM
                % Apply scaling factor to the data
                dataToFit   = op_takesubspec(MRSCont.processed.metab{kk},'A');
                basisSetOff = MRSCont.fit.basisSet;
                basisSetOff.names{1} = 'A'; 
                basisSetOff.fids = basisSetOff.fids(:,:,1);
                basisSetOff.specs = basisSetOff.specs(:,:,1);
                dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
                fitOpts.GAP = MRSCont.opts.fit.GAP.A;
    
    
                % Call the fit function
                [fitParamsOff, resBasisSetOff] = fit_runFit(dataToFit, basisSetOff, fitModel, fitOpts);
    
                % Save back the basis set and fit parameters to MRSCont
                MRSCont.fit.resBasisSet.metab{1,kk,1}             = resBasisSetOff;
                MRSCont.fit.results.metab.fitParams{1,kk,1}   = fitParamsOff;
    
    
                %%% 2b. FIT DIFF1-SPECTRUM
                % Apply scaling factor to the data
                fitOpts     = MRSCont.opts.fit;
                dataToFit   = op_takesubspec(MRSCont.processed.metab{kk},'diff1');
                basisSetDiff1 = MRSCont.fit.basisSet;
                basisSetDiff1.fids = basisSetDiff1.fids(:,:,3);
                basisSetDiff1.specs = basisSetDiff1.specs(:,:,3);
                dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
                dataToFit.refShift   = fitParamsOff.refShift;
                dataToFit.refFWHM   = fitParamsOff.refFWHM;
    
                if ~strcmp(fitOpts.coMM3, 'none')
                    fitOpts.CrFactor = 1;
                    [basisSetDiff1] = osp_addDiffMMPeaks(basisSetDiff1,basisSetOff,fitOpts);
                end
    
                fitOpts.GAP = MRSCont.opts.fit.GAP.diff1;
                basisSetDiff1.names{1} = 'diff1';
                
                
                % Automatically pick the metabolties to be included based on
                % the editing target metabolite is set to 100% and only
                % metabolites > 20% co-editing are included

                intensity = sum(abs(basisSetDiff1.specs(basisSetDiff1.ppm > fitOpts.range(1) & basisSetDiff1.ppm < fitOpts.range(2),:)),1);
                idx = find(strcmp(MRSCont.opts.editTarget{1},basisSetDiff1.name));
                if ~isempty(idx)
                    rel_intensity = intensity/intensity(idx);
                    editing_thresh = find(rel_intensity > 0.2);
                    thresh_MetabList = basisSetDiff1.name(editing_thresh);
                    metabList = fit_createMetabList(thresh_MetabList);
                    %Create the modified basis set
                    basisSetDiff1 = fit_selectMetabs(basisSetDiff1, metabList, 1);
                else
                    if strcmp(MRSCont.opts.editTarget{1},'GABA')
                        metabList = fit_createMetabList({'GABA','GSH','Gln','Glu','NAAG','NAA','MM09'});
                    end
                end

                % For manual overwrite indicate the metabolites to be added
                % metabList = fit_createMetabList({'GABA','GSH','Gln','Glu','NAAG','NAA','MM09'});
                % basisSetDiff1 = fit_selectMetabs(basisSetDiff1, metabList, 1);

                

                % Call the fit function
                [fitParamsDiff1, resBasisSetDiff1]  = fit_runFit(dataToFit, basisSetDiff1, fitModel, fitOpts);
    
                % Save back the basis set and fit parameters to MRSCont
                MRSCont.fit.resBasisSet.metab{1,kk,2}           = resBasisSetDiff1;
                MRSCont.fit.results.metab.fitParams{1,kk,2} = fitParamsDiff1;
    
                %Modeling MM spectra after the main spectrum
                if MRSCont.flags.hasMM == 1
                    dataToFit_mm   = op_takesubspec(MRSCont.processed.mm{kk},'A');
                    dataToFit_mm   = op_ampScale(dataToFit_mm, 1/MRSCont.fit.scale{kk});
                     [refShift_mm, ~] = fit_OspreyReferencingMM(dataToFit_mm);
                    [dataToFit_mm]             = op_freqshift(dataToFit_mm,-refShift_mm);
                    %add some info from the metabolite fit
                    dataToFit_mm.lineShape  = fitParamsOff.lineShape;
                    dataToFit_mm.refFWHM  = fitParamsOff.refFWHM;
                    dataToFit_mm.refShift  = 0;
                    % Extract fit options
                    fitOpts_mm    = MRSCont.opts.fit;
                    fitModel_mm    = fitOpts.method;
                    fitOpts_mm.sequence = 'unedited';
    
                     % Call the fit function
                    [fitParams_mm, resBasisSet_mm] = fit_runFitMM(dataToFit_mm, basisSet_mm_A, fitModel_mm, fitOpts_mm);
                     % Save back the basis set and fit parameters to MRSCont
                    MRSCont.fit.basisSet_mm                    = basisSet_mm_A;
                    MRSCont.fit.resBasisSet.mm{1,kk,1}         = resBasisSet_mm;
                    MRSCont.fit.results.mm.fitParams{1,kk,1}   = fitParams_mm;
    
                    % Create full spectralwidth clean MM spectrum
                    disp('Create clean metabolite-nulled  off spectrum...');
                    [mm_clean] = fit_OspreyCleanMM(dataToFit_mm, resBasisSet_mm, fitOpts_mm,MRSCont.fit.scale{kk},fitParams_mm);
                    mm_clean.names = {'A_clean'};
                    mm_clean.flags.isUnEdited = 1;
                    mm_clean.flags.isMEGA = 0;
    %                 [mm_clean,~]          = op_autophase(mm_clean,0.8,1);
                    [refShift_mm, ~] = fit_OspreyReferencingMM(mm_clean);
                    [mm_clean]             = op_freqshift(mm_clean,-refShift_mm);
                    mm_clean   = op_ampScale(mm_clean, MRSCont.fit.scale{kk});
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},mm_clean);
                    dataToFit_mm   = op_takesubspec(MRSCont.processed.mm{kk},'A');
                    [dataToFit_mm]             = op_freqshift(dataToFit_mm,-refShift_mm);
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},dataToFit_mm);
    
    
                    % Create full spectralwidth clean MM spectrum spline model
                    disp('Create clean metabolite-nulled off spline model...');
                    %Fit them with a spline?
                    dataToFit_mm_clean      = op_takesubspec(MRSCont.processed.mm{kk},'A_clean');
                    fitOpts_mm_clean    = MRSCont.opts.fit;
                    fitOpts_mm_clean.range=[dataToFit_mm_clean.ppm(1) dataToFit_mm_clean.ppm(end)];
                    [~,mm_clean_spline] = fit_Osprey_SplineOnly(dataToFit_mm_clean, 0.1, fitOpts_mm_clean.range);
                    mm_clean_spline.names = {'A_spline'};
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},mm_clean_spline);
    
    %                 % Gaussian fit with(out) baseline
    %                 disp('Fit clean metabolite-nulled spectrum...');
    %                 fitOpts_mm.Params2 = fitParams_mm;
    %                 dataToFit_mm_clean   = op_ampScale(dataToFit_mm_clean, 1/MRSCont.fit.scale{kk});
    %                 [fitParams_mm_Gauss, resBasisSet_mm_Gauss] = fit_runFitMM(dataToFit_mm_clean, basisSetMM,  'OspreyMMGaussian', fitOpts_mm);
    %                 MRSCont.fit.basisSet_mm_Gauss                   = basisSetMM;
    %                 MRSCont.fit.resBasisSet.diff1_mm_Gauss{kk}         = resBasisSet_mm_Gauss;
    %                 MRSCont.fit.results.diff1_mm.fitParams{kk}.fitParams_mm_Gauss   = fitParams_mm_Gauss;
    
    %                 % Fit off with individual MM model
                    disp('Fit off with metabolite-nulled spline model...');
                    fitOpts.mm_clean_spline     = mm_clean_spline;
                    fitOpts.prelimParams = fitParamsOff.prelimParams;
                    dataToFit   = op_takesubspec(MRSCont.processed.metab{kk},'A');
                    dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
                    % Call the fit function
                    [fitParamsOffExpMM, resBasisSetOffExpMM]  = fit_runFit(dataToFit, basisSetOff, fitModel, fitOpts);
                    MRSCont.fit.resBasisSet.metab{2,kk,1}         = resBasisSetOffExpMM;
                    MRSCont.fit.results.metab.fitParams{2,kk,1}   = fitParamsOffExpMM;
    
                    %Clean diff1 mm spectrum
                    dataToFit_mm   = op_takesubspec(MRSCont.processed.mm{kk},'diff1');
                    dataToFit_mm   = op_ampScale(dataToFit_mm, 1/MRSCont.fit.scale{kk});
                     [refShift_mm, ~] = fit_OspreyReferencingMM(dataToFit_mm);
                    [dataToFit_mm]             = op_freqshift(dataToFit_mm,-refShift_mm);
                    %add some info from the metabolite fit
                    dataToFit_mm.lineShape  = fitParamsDiff1.lineShape;
                    dataToFit_mm.refFWHM  = fitParamsDiff1.refFWHM;
                    dataToFit_mm.refShift  = 0;
                    % Extract fit options
                    fitOpts_mm    = MRSCont.opts.fit;
                    fitModel_mm    = fitOpts.method;
                    fitOpts_mm.sequence = 'MEGA';
    
                     % Call the fit function
                    [fitParams_mm, resBasisSet_mm] = fit_runFitMM(dataToFit_mm, basisSet_mm_diff1, fitModel_mm, fitOpts_mm);
                     % Save back the basis set and fit parameters to MRSCont
                    MRSCont.fit.basisSet_mm                    = basisSet_mm_diff1;
                    MRSCont.fit.resBasisSet.mm{1,kk,2}         = resBasisSet_mm;
                    MRSCont.fit.results.mm.fitParams{1,kk,2}   = fitParams_mm;
    
                    % Create full spectralwidth clean MM spectrum
                    disp('Create clean metabolite-nulled  diff1 spectrum...');
                    [mm_clean] = fit_OspreyCleanMM(dataToFit_mm, resBasisSet_mm, fitOpts_mm,MRSCont.fit.scale{kk},fitParams_mm);
                    mm_clean.names = {'diff1_clean'};
    %                 [mm_clean,~]          = op_autophase(mm_clean,0.8,1);
                    [refShift_mm, ~] = fit_OspreyReferencingMM(mm_clean);
                    [mm_clean]             = op_freqshift(mm_clean,-refShift_mm);
                    mm_clean   = op_ampScale(mm_clean, MRSCont.fit.scale{kk});
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},mm_clean);
                    dataToFit_mm   = op_takesubspec(MRSCont.processed.mm{kk},'diff1');
                    [dataToFit_mm]             = op_freqshift(dataToFit_mm,-refShift_mm);
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},dataToFit_mm);
    
    
                    % Create full spectralwidth clean MM spectrum spline model
                    disp('Create clean metabolite-nulled diff1 spline model...');
                    %Fit them with a spline?
                    dataToFit_mm_clean      = op_takesubspec(MRSCont.processed.mm{kk},'diff1_clean');
                    fitOpts_mm_clean    = MRSCont.opts.fit;
                    fitOpts_mm_clean.range=[dataToFit_mm_clean.ppm(1) dataToFit_mm_clean.ppm(end)];
                    [~,mm_clean_spline] = fit_Osprey_SplineOnly(dataToFit_mm_clean, 0.1, fitOpts_mm_clean.range);
                    mm_clean_spline.names = {'diff1_spline'};
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},mm_clean_spline);
    
    %                 % Gaussian fit with(out) baseline
    %                 disp('Fit clean metabolite-nulled spectrum...');
    %                 fitOpts_mm.Params2 = fitParams_mm;
    %                 dataToFit_mm_clean   = op_ampScale(dataToFit_mm_clean, 1/MRSCont.fit.scale{kk});
    %                 [fitParams_mm_Gauss, resBasisSet_mm_Gauss] = fit_runFitMM(dataToFit_mm_clean, basisSetMM,  'OspreyMMGaussian', fitOpts_mm);
    %                 MRSCont.fit.basisSet_mm_Gauss                   = basisSetMM;
    %                 MRSCont.fit.resBasisSet.diff1_mm_Gauss{kk}         = resBasisSet_mm_Gauss;
    %                 MRSCont.fit.results.diff1_mm.fitParams{kk}.fitParams_mm_Gauss   = fitParams_mm_Gauss;
    
    %                 % Fit diff1 with individual MM model
                    disp('Fit diff1 with metabolite-nulled spline model...');
                    fitOpts.mm_clean_spline     = mm_clean_spline;
                    fitOpts.prelimParams = fitParamsDiff1.prelimParams;
                    fitOpts.coMM3  = 'none';
                    dataToFit   = op_takesubspec(MRSCont.processed.metab{kk},'diff1');
                    dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
                    % Call the fit function
                    [fitParamsDiff1ExpMM, resBasisSetDiff1ExpMM]  = fit_runFit(dataToFit, basisSetDiff1, fitModel, fitOpts);
                    MRSCont.fit.resBasisSet.metab{2,kk,2}         = resBasisSetDiff1ExpMM;
                    MRSCont.fit.results.metab.fitParams{2,kk,2}   = fitParamsDiff1ExpMM;
    
                end
    
            end
        end
    
    
        %%% 3. CONCATENATED FIT %%%
        % For the concatenated MEGA fit, model the DIFF1 and SUM spectra
        % together.
        if strcmp(fitStyle, 'Concatenated')
            if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.conc.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
                %%% 3a. FIT CONCATENATED SPECTRUM
                % Select the difference and sum spectrum to put into the
                % concatenated fit
                dataToFit   = {op_takesubspec(MRSCont.processed.metab{kk},3) op_takesubspec(MRSCont.processed.metab{kk},4)};
                basisSetConc = MRSCont.fit.basisSet;
                % Create the basis set with difference and sum basis functions
                basisSetConc.fids(:,:,1) = basisSetConc.fids(:,:,3);
                basisSetConc.specs(:,:,1) = basisSetConc.specs(:,:,3);
                basisSetConc.fids(:,:,2) = basisSetConc.fids(:,:,4);
                basisSetConc.specs(:,:,2) = basisSetConc.specs(:,:,4);
                basisSetConc.fids(:,:,3:4) = [];
                basisSetConc.specs(:,:,3:4) = [];
    
                if ~strcmp(fitOpts.coMM3, 'none')
                    basisSetDiff1 = MRSCont.fit.basisSet;
                    basisSetDiff1.fids = basisSetDiff1.fids(:,:,3);
                    basisSetDiff1.specs = basisSetDiff1.specs(:,:,3);
                    fitOpts.CrFactor = 1;
                    [basisSetDiff1] = osp_addDiffMMPeaks(basisSetDiff1,fitOpts);
                    basisSetConc.fids(:,:,1) = basisSetDiff1.fids(:,:);
                    basisSetConc.specs(:,:,1) = basisSetDiff1.specs(:,:);
                end
    
                % Apply scaling factor to the data
                for rr = 1:length(dataToFit)
                    dataToFit{rr}   = op_ampScale(dataToFit{rr}, 1/MRSCont.fit.scale{kk});
                end
                % Call the multi-spectrum fit function
                [fitParamsConc, resBasisSetConc] = fit_runFitMultiplex(dataToFit, basisSetConc, fitModel, fitOpts);
    
                % Save back the basis set and fit parameters to MRSCont
                MRSCont.fit.resBasisSet.conc{kk}           = resBasisSetConc;
                MRSCont.fit.results.conc.fitParams{kk} = fitParamsConc;
            end
        end
        % ----- LCModel wrapper fit pipeline -----
    elseif strcmpi(MRSCont.opts.fit.method, 'LCModel')
         [MRSCont] = LCModelWrapper(MRSCont,kk,progressText);
    end
end
%Fit again with mean MM basis function
if MRSCont.opts.fit.MeanMM
    if MRSCont.flags.hasMM == 1

        fids = zeros(mm_clean_spline.sz(1),MRSCont.nDatasets(1));
        for kk = 1:MRSCont.nDatasets(1)
            fids(:,kk)     = MRSCont.processed.mm{kk}.fids(:,6);
        end
        mean_fids = mean(fids,2);
        mm_clean_spline.fids = mean_fids;
        mm_clean_spline.specs = fftshift(fft(mean_fids,[],1),1);
        fitOpts.mm_clean_spline     = mm_clean_spline;
        MRSCont.opts.fit.mm_clean_spline= mm_clean_spline;
    end

    for kk = 1:MRSCont.nDatasets(1)
        [~] = printLog('OspreyFit',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
         % Fit diff1 with individual MM model
        disp('Fit off with mean metabolite-nulled spline model...');
        % Call the fit function
        fitOpts.prelimParams = MRSCont.fit.results.metab.fitParams{1,kk,1}.prelimParams;
        fitOpts.GAP = MRSCont.opts.fit.GAP.A;

        dataToFit   = op_takesubspec(MRSCont.processed.metab{kk},'A');
        dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
        dataToFit.refShift   = MRSCont.fit.results.metab.fitParams{1,kk,1}.refShift;
        dataToFit.refFWHM   = MRSCont.fit.results.metab.fitParams{1,kk,1}.refFWHM;

        [fitParamsOffExpMM, resBasisSetOffExpMM]  = fit_runFit(dataToFit, basisSetOff, fitModel, fitOpts);
        MRSCont.fit.resBasisSet.metab{3,kk,1}         = resBasisSetOffExpMM;
        MRSCont.fit.results.metab.fitParams{3,kk,1}   = fitParamsOffExpMM;

    end

    if MRSCont.flags.hasMM == 1

        fids = zeros(mm_clean_spline.sz(1),MRSCont.nDatasets(1));
        for kk = 1:MRSCont.nDatasets(1)
            fids(:,kk)     = MRSCont.processed.mm{kk}.fids(:,7);
        end
        mean_fids = mean(fids,2);
        mm_clean_spline.fids = mean_fids;
        mm_clean_spline.specs = fftshift(fft(mean_fids,[],1),1);
        fitOpts.mm_clean_spline     = mm_clean_spline;
        MRSCont.opts.fit.mm_clean_spline= mm_clean_spline;
    else
        load('/Volumes/Samsung/working/editedMM/mm_clean_spline.mat') % This is just for now :D
        [refShift_mm, ~] = fit_OspreyReferencingMM(mm_clean_spline);
        [mm_clean_spline]             = op_freqshift(mm_clean_spline,-refShift_mm);
        MRSCont.opts.fit.mm_clean_spline= mm_clean_spline;
        fitOpts.mm_clean_spline = mm_clean_spline;
    end

    fitOpts.coMM3  = 'none';
    for kk = 1:MRSCont.nDatasets(1)
        [~] = printLog('OspreyFit',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
         % Fit diff1 with individual MM model
        disp('Fit diff1 with mean metabolite-nulled spline model...');
        % Call the fit function
        fitOpts.prelimParams = MRSCont.fit.results.metab.fitParams{1,kk,2}.prelimParams;
        fitOpts.GAP = MRSCont.opts.fit.GAP.diff1;

        dataToFit   = op_takesubspec(MRSCont.processed.metab{kk},'diff1');
        dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
        dataToFit.refShift   = MRSCont.fit.results.metab.fitParams{1,kk,2}.refShift;
        dataToFit.refFWHM   = MRSCont.fit.results.metab.fitParams{1,kk,2}.refFWHM;

        [fitParamsDiff1ExpMM, resBasisSetDiff1ExpMM]  = fit_runFit(dataToFit, basisSetDiff1, fitModel, fitOpts);
        MRSCont.fit.resBasisSet.metab{3,kk,2}         = resBasisSetDiff1ExpMM;
        MRSCont.fit.results.metab.fitParams{3,kk,2}   = fitParamsDiff1ExpMM;

    end
end


time = toc(metFitTime);
[~] = printLog('done',time,MRSCont.nDatasets,1,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
MRSCont.runtime.FitMet = time;

end
