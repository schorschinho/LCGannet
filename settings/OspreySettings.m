function [MRSCont] = OspreySettings
%% [MRSCont] = OspreySettings
%   This function initialises default settings and variables required to
%   run Osprey.
%
%   USAGE:
%       [MRSCont] = OspreySettings;
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-18)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-02-19: First version of the code.

close all;

%%% 1. REQUIRED INPUT FLAGS %%%
% Initialise MRSCont struct with default settings
MRSCont = struct;
MRSCont.flags.isUnEdited            = 0;
MRSCont.flags.isMEGA                = 0;
MRSCont.flags.isHERMES              = 0;
MRSCont.flags.isHERCULES            = 0;
MRSCont.flags.isPRIAM               = 0;
MRSCont.flags.isMRSI                = 0;
MRSCont.flags.isSPECIAL             = 0;
MRSCont.flags.isSERIES              = 0;
MRSCont.flags.addImages             = 0;
MRSCont.flags.reordered             = 0;
MRSCont.opts.savePDF                = 0;
MRSCont.opts.saveLCM                = 0;
MRSCont.opts.savejMRUI              = 0;
MRSCont.opts.saveNII                = 1;
MRSCont.opts.saveVendor             = 0;
MRSCont.opts.fit.includeMetabs      = {'default'};      % Options: 'default', 'full', custom set
MRSCont.opts.fit.method             = 'Osprey';         % Options: 'Osprey' (default), 'LCModel'
MRSCont.opts.fit.range              = [0.2 4.2];        % Default: [0.2 4.2]
MRSCont.opts.fit.rangeWater         = [2.0 7.4];        % Default: [2.0 7.4]
MRSCont.opts.fit.style              = 'Separate';       % Options: 'Separate' (Default - will fit DIFF and OFF separately), 'Concatenated' (will fit DIFF and SUM simultaneously),
MRSCont.opts.fit.bLineKnotSpace     = 0.4;              % Baseline spline knot spacing [ppm]. Default: 0.4.
MRSCont.opts.fit.fitMM              = 1;                % Add MM and lipid basis functions to basis set? Default: 1.
MRSCont.opts.fit.coMM3              = 'none';           % Add co-edited MM3 peak model for GABA editing? Default: none.
MRSCont.opts.fit.FWHMcoMM3          = 14;               % FWHM [Hz] of the co-edited peak Default: 14 Hz.
MRSCont.opts.fit.RegParLS           = 1000;             % Regularization parameter of arbitrary lineshape conv kernel
MRSCont.opts.ECC.raw                = 1;                % Do ECC for all metabolite spectra.
MRSCont.opts.ECC.mm                 = 1;                 % Do ECC for all metabolite-nulled spectra.
MRSCont.opts.cosmetics.LB           = 0;                % Do cosmetic LB
MRSCont.opts.cosmetics.Zoom         = 2.75;             % Do cosmetic Zoom
MRSCont.opts.img.deface             = 0;                % Deface data
MRSCont.opts.load.undoPhaseCycle    = 1;                % Undo phase cycle for nifti MRS
%%% 2. FIND AND SET PATHS %%%
% Osprey
[settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
if isempty(settingsFolder)
     error('Osprey not found! Please install Osprey (https://github.com/schorschinho/osprey) and include it in your MATLAB path.');
else
    allFolders      = strsplit(settingsFolder, filesep);
    ospFolder       = strjoin(allFolders(1:end-1), filesep); % parent folder (= Osprey folder)
end

% SPM
% Check if SPM12 is installed
spmversion = fileparts(which('spm'));
if isempty(spmversion)
    warning('SPM not found! Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and include it in your MATLAB path.');
elseif strcmpi(spmversion(end-3:end),'spm8')
    warning(['SPM8 detected, but only SPM12 is supported. ' ...
           'Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and include it in your MATLAB path.']);
end

%%% 3. INITIALISE MRSCONT FIELDS AND FLAGS %%%
% Set default fields
MRSCont.ospFolder           = ospFolder;
MRSCont.files               = {};
MRSCont.files_mm            = {};
MRSCont.files_ref           = {};
MRSCont.files_w             = {};
MRSCont.files_mm_ref        = {};
MRSCont.flags.hasMM         = 0;
MRSCont.flags.hasRef        = 0;
MRSCont.flags.hasWater      = 0;
MRSCont.flags.hasMMRef      = 0;
% Set default flags
MRSCont.flags.didLCMWrite   = 0;
MRSCont.flags.didjMRUIWrite = 0;
MRSCont.flags.didVendorWrite= 0;
MRSCont.flags.didJob        = 0;
MRSCont.flags.didLoad   = 0;
MRSCont.flags.didProcess    = 0;
MRSCont.flags.didCoreg      = 0;
MRSCont.flags.didSeg        = 0;
MRSCont.flags.didFit        = 0;
MRSCont.flags.didQuantify   = 0;
MRSCont.flags.didOverview   = 0;

if (ismcc || isdeployed)
    if ismac
        currentDir = ctfroot;
        [currentDir,~,~] = fileparts(currentDir);
         SepFileList =  split(currentDir, filesep);
        index = find(strcmp(SepFileList,'application'));
        if ~isempty(index)
            MRSCont.opts.fit.basissetFolder = fullfile(SepFileList{1:index},'basissets');
        else
            MRSCont.opts.fit.basissetFolder = [];
        end
    end
    if isunix
    	currentDir = pwd;
    	SepFileList =  split(currentDir, filesep);
        index = find(strcmp(SepFileList,'application'));
        if ~isempty(index)
            MRSCont.opts.fit.basissetFolder = fullfile(SepFileList{1:index},'basissets');
        else
            MRSCont.opts.fit.basissetFolder = [];
        end
    end
    if ispc
        [~, result] = system('path');
        currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
        SepFileList =  split(currentDir, filesep);
        MRSCont.opts.fit.basissetFolder =currentDir;
        index = find(strcmp(SepFileList,'application'));
        if ~isempty(index)
            MRSCont.opts.fit.basissetFolder = fullfile(SepFileList{1:index},'basissets');
        else
            MRSCont.opts.fit.basissetFolder = [];
        end
    end
    
else
    MRSCont.opts.fit.basissetFolder = [];
end


end
