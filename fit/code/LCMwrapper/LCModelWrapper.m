function [MRSCont] = LCModelWrapper(MRSCont,kk,progressText)
%% [MRSCont] = LCModelWrapper(MRSCont)
%   This function performs spectral fitting in LCModel.
%
%   USAGE:
%       [MRSCont] = osp_fitUnEdited(MRSCont, basisSet);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       basisSet    = Osprey basis set container.
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
%       2023-07-12: First version of the code.

    [~] = printLog('OspreyFit',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);

    % Put the scale to be 1 for now. Might have to adjust.
    MRSCont.fit.scale{kk} = 1;

    % Create the sub-folder where the LCModel results will be saved,
    % otherwise LCModel will throw 'FATAL ERROR MAIN 2'.
    if ~exist(fullfile(MRSCont.outputFolder, 'LCMoutput'), 'dir')
        mkdir(fullfile(MRSCont.outputFolder, 'LCMoutput'));
    end


    % Find the path to LCModel binary if no custom path has been provided
    if ~isfield(MRSCont.opts.fit, 'customLCModelBinary') 

        % Assemble path according to system architecture
        pathLCModelBinaryStruct = osp_platform('lcmodel');
        pathLCModelBinary = fullfile(MRSCont.ospFolder, 'libraries', 'LCModel',pathLCModelBinaryStruct.os,pathLCModelBinaryStruct.osver,['LCModel_' pathLCModelBinaryStruct.os,'_',pathLCModelBinaryStruct.osver]);
        switch  pathLCModelBinaryStruct.os
            case 'macos'
                bin = 'lcmodel';
            case 'unix'
                bin = 'lcmodel';
            case 'win'
                bin = 'LCModel.exe';
        end

        pathLCModelBinaryFile = fullfile(pathLCModelBinary, bin);

        if ~exist(pathLCModelBinaryFile,'file') %Unzip binary
            unzip([pathLCModelBinaryFile '.zip'], [pathLCModelBinary filesep]);
        end

        if ~exist(pathLCModelBinaryFile,'file') %Binary is still missing
            error('ERROR: No LCModel binary found.  Aborting!!');
        else
            if ~(ismcc || isdeployed)
                addpath(which('libraries/LCModel'));
            end
        end

    else

        % If custom path has been provided, check that it is valid and the binary
        % exists
        pathLCModelBinaryFile = MRSCont.opts.fit.customLCModelBinary;
        if ~exist(pathLCModelBinaryFile,'file') %Binary is still missing
            error('ERROR: Custom LCModel binary path is not a valid file.  Aborting!!');
        else
            if ~(ismcc || isdeployed)
                addpath(which('libraries/LCModel'));
            end
        end

    end

    % check that the binary can in fact be executed
    [~,attribs] = fileattrib(pathLCModelBinaryFile);

    if (~isnan(attribs.UserExecute) && ~attribs.UserExecute)
        disp('LCModel binary is not executable; attempting to correct')
        try
            fileattrib(pathLCModelBinaryFile,'+x','u')
            % this might fail on certain platforms, certain filesystems, or with unexpected ownership
        catch
            disp('...sadly, that did not work.')
        end
    end

    % Call LCModel and read in the LCModel output files
    if MRSCont.flags.isUnEdited
        callLCModel(MRSCont.opts.fit.lcmodel.controlfileA{kk}, pathLCModelBinaryFile);
        MRSCont.fit.results.metab.fitParams{1,kk,1} = readLCMFitParams(MRSCont, 'A', kk);
    elseif MRSCont.flags.isMEGA
        callLCModel(MRSCont.opts.fit.lcmodel.controlfileA{kk}, pathLCModelBinaryFile);
        MRSCont.fit.results.metab.fitParams{1,kk,1} = readLCMFitParams(MRSCont, 'A', kk);
        callLCModel(MRSCont.opts.fit.lcmodel.controlfileDiff1{kk}, pathLCModelBinaryFile);
        MRSCont.fit.results.metab.fitParams{1,kk,2} = readLCMFitParams(MRSCont, 'diff1', kk);
    end
end
