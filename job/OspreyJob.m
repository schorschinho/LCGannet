function [MRSCont] = OspreyJob(jobFile,GUI,debug)
%% [MRSCont] = OspreyJob(jobFile)
%   This function loads the Osprey job defined in jobFile.
%   The job can be submitted in the following formats:
%       - .m
%
%   A valid Osprey job contains four distinct classes of items:
%       1. basic information on the MRS sequence used
%       2. several settings for data handling and modeling
%       3. a list of MRS (and, optionally, structural imaging) data files
%          to be loaded
%       4. an output folder to store the results and exported files
%
%   See the example files in the /exampledata/ folder for details.
%
%   USAGE:
%       [MRSCont] = OspreyJob(jobFile);
%
%   INPUTS:
%       jobFile     = File containing a correct Osprey job definition.
%                     Accepted file formats are .m and .csv.
%       GUI         = flag to decide whether plot is used in GUI
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-07-15)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-07-15: First version of the code.


% Close any remaining open figures and parse input arguments
close all;
warning('off','all');
if nargin < 3
    debug = '00';
    if nargin<2
        GUI = 0;
    end
end

%%% 1. INITIALISE DATA CONTAINER WITH DEFAULT SETTINGS
[MRSCont] = OspreySettings;
VersionStruct  = getCurrentVersion;

%%% 2. CHECK JOB INPUT FILE FORMAT %%%
[~,~,ext] = fileparts(jobFile);
switch ext
    case '.m'
        jobFileFormat = 'm';
    case '.csv'
        jobFileFormat = 'csv';
    case '.json'
        jobFileFormat = 'json';
    otherwise
        error('Unrecognized job file datatype. Job files need to end in .CSV, .M or .JSON');
end


%%% 3a. LOAD AND RUN M FILE %%%
if strcmp(jobFileFormat, 'm')
    run(jobFile);
end

%%% 3b. LOAD CSV FILE %%%
if strcmp(jobFileFormat,'csv')
    jobCSV = readtable(jobFile, 'Delimiter', ',');
    % UTF-8 encoding may screw up the first column name; correct here
    if ~strcmp(jobCSV.Properties.VariableNames{1}, 'files')
        jobCSV.Properties.VariableNames{1} = 'files';
    end
    % Convert to a struct array
    jobStruct = table2struct(jobCSV);

    % Check whether the relevant fieldnames have been entered,
    % and save them as separate cells to be saved into the MRSCont
    % container.
    if isfield(jobStruct, 'files')
        files = {jobStruct.files};
    else
        error('Invalid job file! A job file needs to contain at least metabolite data in the field ''files''.');
    end
    if isfield(jobStruct, 'files_mm')  %re_mm Adding functionality for MM
        files_mm = {jobStruct.files_mm};   %re_mm
    end %re_mm
    if isfield(jobStruct, 'files_mm_ref')
        files_mm_ref = {jobStruct.files_mm_ref};
    end
    if isfield(jobStruct, 'files_ref')
        files_ref = {jobStruct.files_ref};
    end
    if isfield(jobStruct, 'files_w')
        files_w = {jobStruct.files_w};
    end
    if isfield(jobStruct, 'files_nii')
        files_nii = {jobStruct.files_nii};
    end
    if isfield(jobStruct, 'files_seg')
        files_seg = {jobStruct.files_seg};
    end
    if isfield(jobStruct, 'files_sense')
        files_sense = {jobStruct.sense};
    end
    if isfield(jobStruct, 'outputFolder')
        outputFolder = jobStruct(1).outputFolder;
    else
        error('Invalid job file! A job file needs to specify an output folder.');
    end
    if isfield(jobStruct,'seqType')
        seqType = jobStruct(1).seqType;
    else
        fprintf('Sequence type is set to unedited (default). Please indicate otherwise in the csv-file or the GUI \n');
        seqType = 'unedited';
    end
    if isfield(jobStruct,'MultiVoxel')
        MultiVoxel = jobStruct(1).MultiVoxel;
    else
        fprintf('Multivoxel is set to single voxel spectroscopy  (default). Please indicate otherwise in the csv-file or the GUI \n');
        MultiVoxel = 'SVS';
    end
    if isfield(jobStruct,'editTarget')
        MRSCont.opts.editTarget = jobStruct(1).editTarget;
    else
        fprintf('Editing target is set to none (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.editTarget = 'unedited';
    end
    if isfield(jobStruct,'savePDF')
        MRSCont.opts.savePDF = jobStruct(1).savePDF;
    else
        fprintf('PDF-output will not be stored (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.savePDF = 0;
    end
    if isfield(jobStruct,'dataScenario')
        dataScenario = jobStruct(1).dataScenario;
    else
        fprintf('Data scenario is set to ''invivo'' (default). Please indicate otherwise in the csv-file or the GUI \n');
        dataScenario = 'invivo';
    end
    if isfield(jobStruct,'SpecReg')
        MRSCont.opts.SpecReg = jobStruct(1).SpecReg;
    else
        fprintf('Spectral Registration is set to ''RobSpecReg'' (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.SpecReg = 'RobSpecReg';
    end
    if isfield(jobStruct,'saveLCM')
        MRSCont.opts.saveLCM = jobStruct(1).saveLCM;
    else
        fprintf('LCModel-readable files will be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.saveLCM = 1;
    end
    if isfield(jobStruct,'savejMRUI')
        MRSCont.opts.savejMRUI = jobStruct(1).savejMRUI;
    else
        fprintf('jMRUI-readable files will be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.savejMRUI = 1;
    end
    if isfield(jobStruct,'saveVendor')
        MRSCont.opts.saveVendor = jobStruct(1).saveVendor;
    else
        fprintf('Vendor-specific files (SDAT/SPAR, RDA, P) will be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.saveVendor = 1;
    end
    if isfield(jobStruct,'saveNII')
        MRSCont.opts.saveNII = jobStruct(1).saveNII;
    else
        fprintf('NIfTI-MRS files will not be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.saveNII = 0;
    end
    if isfield(jobStruct,'includeMetabs')
        opts.fit.includeMetabs = jobStruct(1).includeMetabs;
    else
        fprintf('Included metabolites is set to default (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.fit.includeMetabs = {'default'};
    end
    if isfield(jobStruct,'method')
        MRSCont.opts.fit.method = jobStruct(1).method;
    else
        fprintf('Fitting algorithm is set to Osprey (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.fit.method = 'Osprey';
    end
    if isfield(jobStruct,'style')
        opts.fit.style = jobStruct(1).style;
    else
        fprintf('Fitting style is set to Concatenated (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.fit.style = 'Concatenated';
    end
    if isfield(jobStruct,'lolim_range') && isfield(jobStruct,'uplim_range')
        opts.fit.range = [jobStruct(1).lolim_range jobStruct(1).uplim_range];
    else
        fprintf('Fitting range is set to [0.5 4] ppm (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.fit.range = [0.5 4] ;
    end
    if isfield(jobStruct,'lolim_rangew') && isfield(jobStruct,'uplim_rangew')
        MRSCont.opts.fit.rangeWater = [jobStruct(1).lolim_rangew jobStruct(1).uplim_rangew];
    else
        fprintf('Fitting range is set to [2.0 7.4] ppm (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.fit.rangeWater = [2.0 7.4] ;
    end
    if isfield(jobStruct,'KnotSpace')
        MRSCont.opts.fit.bLineKnotSpace = jobStruct(1).KnotSpace;
    else
        fprintf('Baseline knot spacing is set to 0.4 ppm (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.fit.bLineKnotSpace = 0.4;
    end
    if isfield(jobStruct,'fitMM')
        MRSCont.opts.fit.fitMM = jobStruct(1).fitMM;
    else
        fprintf('Adding macromolecule and lipid basis functions to the fit (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.fit.fitMM = 1;
    end
    if isfield(jobStruct,'RelaxationAtlas')
        MRSCont.opts.quantify.RelaxationAtlas = jobStruct(1).RelaxationAtlas;
    else
        MRSCont.opts.quantify.RelaxationAtlas = 0;
    end
    if isfield(jobStruct,'RelaxationAtlasAge')
        MRSCont.opts.quantify.RelaxationAtlasAge = jobStruct(1).RelaxationAtlasAge;
    else
        MRSCont.opts.quantify.RelaxationAtlasAge = 0;
    end
    if isfield(jobStruct,'deface')
        MRSCont.opts.img.deface = jobStruct.deface;
    else
        fprintf('Structrual images are not defaced (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.img.deface = 0;
    end
    if isfield(jobStruct,'exportParams')
        MRSCont.opts.exportParams.flag = 1;
        MRSCont.opts.exportParams.path = jobStruct(1).exportParams;
        assert(isstring(jobStruct(1).exportParams), 'When exporting spectral fitting parameters, "exportParams" must be a string specifying the save directory')
    else
        fprintf('Spectral fitting parameters will not be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.exportParams.flag = 0;
        MRSCont.opts.exportParams.path = {};
    end
end

%%% 3b. LOAD JSON FILE %%%
if strcmp(jobFileFormat,'json')
    fid = fopen(jobFile); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    if strcmp('win',osp_platform('filesys'))
        str = strrep(str,'\','\\');
    end
    pattern = '[ \t\n]*“'; % Match zero or more spaces, tabs, or newlines, followed by a double quote
    replacement = '"'; % Replace the matched string with just a double quote
    str = regexprep(str, pattern, replacement);
    jobStruct  = jsondecode(str);

    % Check whether the relevant fieldnames have been entered,
    % and save them as separate cells to be saved into the MRSCont
    % container.
    if isfield(jobStruct, 'files')
        if ~iscell(jobStruct.files{1}) || size(jobStruct.files,1) == 1
            files = jobStruct.files';
        end
        if iscell(files{1})
            files = cellfun(@transpose,files,'UniformOutput',false);
        end
    else
        error('Invalid job file! A job file needs to contain at least metabolite data in the field ''files''.');
    end
    if isfield(jobStruct, 'files_mm')  %re_mm Adding functionality for MM
        if ~iscell(jobStruct.files_mm{1}) || size(jobStruct.files_mm,1) == 1
            files_mm = jobStruct.files_mm';   %re_mm
        end
        if iscell(files_mm{1})
            files_mm = cellfun(@transpose,files_mm,'UniformOutput',false);
        end
    end %re_mm
    if isfield(jobStruct, 'files_mm_ref')
        if ~iscell(jobStruct.files_mm_ref{1}) || size(jobStruct.files_mm_ref,1) == 1
            files_mm_ref = jobStruct.files_mm_ref';
        end
        if iscell(files_mm_ref{1})
            files_mm_ref = cellfun(@transpose,files_mm_ref,'UniformOutput',false);
        end
    end
    if isfield(jobStruct, 'files_ref')
        if ~iscell(jobStruct.files_ref{1}) || size(jobStruct.files_ref,1) == 1
            files_ref = jobStruct.files_ref';
        end
        if iscell(files_ref{1})
            files_ref = cellfun(@transpose,files_ref,'UniformOutput',false);
        end
    end
    if isfield(jobStruct, 'files_w')
        if ~iscell(jobStruct.files_w{1}) || size(jobStruct.files_w,1) == 1
            files_w = jobStruct.files_w';
        end
        if iscell(files_w{1})
            files_w = cellfun(@transpose,files_w,'UniformOutput',false);
        end
    end
    if isfield(jobStruct, 'files_nii')
        files_nii = jobStruct.files_nii';
    end
    if isfield(jobStruct, 'files_seg')
        files_seg = jobStruct.files_seg';
    end
    if isfield(jobStruct, 'files_sense')
        files_sense = jobStruct.sense';
    end
    if isfield(jobStruct, 'outputFolder')
        outputFolder = jobStruct.outputFolder{1};
    else
        error('Invalid job file! A job file needs to specify an output folder.');
    end
    if isfield(jobStruct, 'file_stat')
        file_stat = jobStruct.file_stat{1};
    else
        '';
    end
    if isfield(jobStruct,'seqType')
        seqType = jobStruct.seqType;
    else
        seqType = 'unedited';
    end
    if isfield(jobStruct,'MultiVoxel')
        MultiVoxel = jobStruct.MultiVoxel;
    else
        MultiVoxel = 'SVS';
    end
    if isfield(jobStruct,'editTarget')
        editTarget = jobStruct.editTarget';
    else
        editTarget = {'none'};
    end
    if isfield(jobStruct,'savePDF')
        MRSCont.opts.savePDF = jobStruct.savePDF;
    else
        MRSCont.opts.savePDF = 0;
    end
    if isfield(jobStruct,'dataScenario')
        dataScenario = jobStruct.dataScenario;
    else
        dataScenario = 'invivo';
    end
    if isfield(jobStruct,'SpecReg')
        MRSCont.opts.SpecReg = jobStruct.SpecReg;
    else
        MRSCont.opts.SpecReg = 'RobSpecReg';
    end
    if isfield(jobStruct,'SubSpecAlignment')
        MRSCont.opts.SubSpecAlignment = jobStruct.SubSpecAlignment;
    else
        MRSCont.opts.SubSpecAlignment = 'L2Norm';
    end
    if isfield(jobStruct,'UnstableWater')
        MRSCont.opts.UnstableWater = str2num(jobStruct.UnstableWater);
    else
        MRSCont.opts.UnstableWater = 0;
    end
    if isfield(jobStruct,'ECCmetab')
        MRSCont.opts.ECC.raw = str2num(jobStruct.ECCmetab);
    end
    if isfield(jobStruct,'ECCmm')
        MRSCont.opts.ECC.mm = str2num(jobStruct.ECCmm);
    else
        MRSCont.opts.ECC.mm = 1;
    end
    if isfield(jobStruct,'MM3coModel')
        opts.fit.coMM3 = jobStruct.MM3coModel;
    end
    if isfield(jobStruct,'FWHMMM3co')
        opts.fit.FWHMcoMM3= str2num(jobStruct.FWHMMM3co);
        if isempty(opts.fit.FWHMcoMM3) || strcmp(opts.fit.FWHMcoMM3,' ')
            opts.fit.FWHMcoMM3 = 14;
        end
    end
    if isfield(jobStruct,'saveLCM')
        MRSCont.opts.saveLCM = str2num(jobStruct.saveLCM);
    else
        MRSCont.opts.saveLCM = 0;
    end
    if isfield(jobStruct,'savejMRUI')
        MRSCont.opts.savejMRUI = str2num(jobStruct.savejMRUI);
    else
        MRSCont.opts.savejMRUI = 0;
    end
    if isfield(jobStruct,'saveVendor')
        MRSCont.opts.saveVendor = str2num(jobStruct.saveVendor);
    else
        MRSCont.opts.saveVendor = 0;
    end
    if isfield(jobStruct,'saveNII')
        MRSCont.opts.saveNII = str2num(jobStruct.saveNII);
    else
        MRSCont.opts.saveNII = 0;
    end
    if isfield(jobStruct,'savePDF')
        MRSCont.opts.savePDF = str2num(jobStruct.savePDF);
    else
        MRSCont.opts.savePDF = 0;
    end
    if isfield(jobStruct,'includeMetabs')
        opts.fit.includeMetabs = jobStruct.includeMetabs';
    else
        opts.fit.includeMetabs = {'default'};
    end
    if isfield(jobStruct,'method')
        opts.fit.method = jobStruct.method;
    else
        opts.fit.method = 'Osprey';
    end
    if isfield(jobStruct,'style')
        opts.fit.style = jobStruct.style;
    else
        opts.fit.style = 'Separate';
    end
    if isfield(jobStruct,'lolim_range') && isfield(jobStruct,'uplim_range')
        opts.fit.range = [str2num(jobStruct.lolim_range) str2num(jobStruct.uplim_range)];
    else
        MRSCont.opts.fit.range = [0.2 4.2] ;
    end
    if isfield(jobStruct,'lolim_rangew') && isfield(jobStruct,'uplim_rangew')
        MRSCont.opts.fit.rangeWater = [str2num(jobStruct.lolim_rangew) str2num(jobStruct.uplim_rangew)];
    else
        MRSCont.opts.fit.rangeWater = [2.0 7.4] ;
    end
    if isfield(jobStruct,'bLineKnotSpace')
        MRSCont.opts.fit.bLineKnotSpace = str2num(jobStruct.bLineKnotSpace);
    else
        MRSCont.opts.fit.bLineKnotSpace = 0.4;
    end
    if isfield(jobStruct,'fitMM')
        MRSCont.opts.fit.fitMM = str2num(jobStruct.fitMM);
    else
        MRSCont.opts.fit.fitMM = 1;
    end
    if isfield(jobStruct,'RelaxationAtlas')
        MRSCont.opts.quantify.RelaxationAtlas = jobStruct.RelaxationAtlas;
    else
        MRSCont.opts.quantify.RelaxationAtlas = 0;
    end
    if isfield(jobStruct,'RelaxationAtlasAge')
        MRSCont.opts.quantify.RelaxationAtlasAge = jobStruct.RelaxationAtlasAge;
    else
        MRSCont.opts.quantify.RelaxationAtlasAge = 0;
    end
    if isfield(jobStruct,'basisSet')
        if isfolder(jobStruct.basisSet)
            opts.fit.basissetFolder = jobStruct.basisSet;
        else
            opts.fit.basisSetFile = jobStruct.basisSet;
        end
    end
    if isfield(jobStruct,'deface')
        opts.img.deface = str2num(jobStruct.deface);
    else
        opts.img.deface = 0;
    end
    if isfield(jobStruct,'undoPhaseCycle')
        opts.load.undoPhaseCycle = str2num(jobStruct.undoPhaseCycle);
    else
        opts.load.undoPhaseCycle = 1;
    end
    if isfield(jobStruct,'SubSpecOrder')
        MRSCont.opts.Order = jobStruct.SubSpecOrder';
    else
        MRSCont.opts.Order = [];
    end
    debug = '11';
    if isfield(jobStruct,'exportParams')
        MRSCont.opts.exportParams.flag = 1;
        MRSCont.opts.exportParams.path = jobStruct.exportParams;
        assert(isstring(jobStruct.exportParams), 'When exporting spectral fitting parameters, "exportParams" must be a string specifying the save directory')
    else
        fprintf('Spectral fitting parameters will not be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        MRSCont.opts.exportParams.flag = 0;
        MRSCont.opts.exportParams.path = {};
    end
end

if exist('opts','var')
    names = fields(opts);
    for f = 1 : length(names)
        if ~strcmp(names{f},'fit') && ~strcmp(names{f},'ECC')
            MRSCont.opts.(names{f}) = opts.(names{f});
        else
            names_second = fields(opts.(names{f}));
            for nf = 1 : length(names_second)
                MRSCont.opts.(names{f}).(names_second{nf}) = opts.(names{f}).(names_second{nf}); 
            end
        end
    end
end

if ~isfield(MRSCont.opts, 'UnstableWater')
    MRSCont.opts.UnstableWater = 0;
end

if ~isfield(MRSCont.opts, 'Order')
    MRSCont.opts.Order = [];
end

if ~isfield(MRSCont.opts, 'SubSpecAlignment')
    MRSCont.opts.SubSpecAlignment.mets = 'L2Norm';
    MRSCont.opts.SubSpecAlignment.mm = 'none';
else if ~isfield(MRSCont.opts.SubSpecAlignment, 'mets')
    SubSpecAlignment = MRSCont.opts.SubSpecAlignment;
    MRSCont.opts = rmfield(MRSCont.opts,'SubSpecAlignment');
    MRSCont.opts.SubSpecAlignment.mets = SubSpecAlignment;
    MRSCont.opts.SubSpecAlignment.mm = 'none'; 
    end
end

if ~isfield(MRSCont.opts,'exportParams')
    fprintf('Spectral fitting parameters will not be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
    MRSCont.opts.exportParams.flag = 0;
    MRSCont.opts.exportParams.path = {};
end


%%% 4. SAVE SETTINGS & STAT FILE INTO MRSCONT  %%%
% Parse the sequence type entry
switch seqType
    case 'unedited'
        MRSCont.flags.isUnEdited    = 1;
        MRSCont.opts.editTarget             = {'none'};
        MRSCont.opts.fit.style = opts.fit.style; 
        if strcmp(opts.fit.style, 'Concatenated')
            fprintf('Fitting style was changed to Separate, because this is unedited data.\n');
            MRSCont.opts.fit.style = 'Separate';
        end 
        if ~isfield(MRSCont.opts.fit, 'GAP')
            MRSCont.opts.fit.GAP.A = [];
            MRSCont.opts.fit.GAP.mm = [];
            MRSCont.opts.fit.GAP.ref = [];
            MRSCont.opts.fit.GAP.ref_mm = [];
            MRSCont.opts.fit.GAP.w = [];
        end
        if ~isfield(MRSCont.opts.fit, 'MeanMM')
            MRSCont.opts.fit.MeanMM = 0;
        end
    case 'MEGA'
        MRSCont.flags.isMEGA        = 1;
        if exist('editTarget','var')
            MRSCont.opts.editTarget             = editTarget;
        else
            MRSCont.opts.editTarget = {'GABA'};
        end
        MRSCont.opts.fit.style = opts.fit.style;
        if strcmp(opts.fit.style, 'Concatenated')
            fprintf('Fitting style was changed to Separate, because concatenated modeling is still under development.\n');
            MRSCont.opts.fit.style = 'Separate';
        end
        if isfield(opts.fit, 'coMM3')
            MRSCont.opts.fit.coMM3 = opts.fit.coMM3;
            MRSCont.opts.fit.FWHMcoMM3 = opts.fit.FWHMcoMM3;
        else
            MRSCont.opts.fit.coMM3 = 'none';
            MRSCont.opts.fit.FWHMcoMM3 = 14;
        end
        if isfield(MRSCont.opts.fit, 'FWHMcoMM3')
            if isempty(MRSCont.opts.fit.FWHMcoMM3)
                MRSCont.opts.fit.FWHMcoMM3 = 14;
            end
        end
        if ~isfield(MRSCont.opts.fit, 'GAP')
            MRSCont.opts.fit.GAP.A = [];
            MRSCont.opts.fit.GAP.diff1 = [];
            MRSCont.opts.fit.GAP.A_mm = [];
            MRSCont.opts.fit.GAP.diff1_mm = [];
            MRSCont.opts.fit.GAP.ref = [];
            MRSCont.opts.fit.GAP.ref_mm = [];
            MRSCont.opts.fit.GAP.w = [];
        else if ~isfield(MRSCont.opts.fit.GAP, 'A')
                MRSCont.opts.fit.GAP.A = [];
            end
            if ~isfield(MRSCont.opts.fit.GAP, 'diff1')
                MRSCont.opts.fit.GAP.diff1 = [];
            end
        end
        if ~isfield(MRSCont.opts.fit, 'MeanMM')
            MRSCont.opts.fit.MeanMM = 0;
        end
    case 'HERMES'
        MRSCont.flags.isHERMES      = 1;
        if exist('editTarget','var')
            MRSCont.opts.editTarget             = editTarget;
        else
            MRSCont.opts.editTarget = {'GABA','GSH'};
        end
        MRSCont.opts.fit.style = opts.fit.style;
        if strcmp(opts.fit.style, 'Concatenated')
            fprintf('Fitting style was changed to Separate, because concatenated modeling is still under development.\n');
            MRSCont.opts.fit.style = 'Separate';
        end
        if isfield(opts.fit, 'coMM3')
            MRSCont.opts.fit.coMM3 = opts.fit.coMM3;
            MRSCont.opts.fit.FWHMcoMM3 = opts.fit.FWHMcoMM3;
        else
            MRSCont.opts.fit.coMM3 = 'fixedGauss';
            MRSCont.opts.fit.FWHMcoMM3 = 14;
        end
        if ~isfield(MRSCont.opts.fit, 'GAP')
            MRSCont.opts.fit.GAP.sum = [];
            MRSCont.opts.fit.GAP.diff1 = [];
            MRSCont.opts.fit.GAP.diff2 = [];
            MRSCont.opts.fit.GAP.diff3 = [];
        else if ~isfield(MRSCont.opts.fit.GAP, 'sum')
                MRSCont.opts.fit.GAP.sum = [];
            end
            if ~isfield(MRSCont.opts.fit.GAP, 'diff1')
                MRSCont.opts.fit.GAP.diff1 = [];
            end
            if ~isfield(MRSCont.opts.fit.GAP, 'diff2')
                MRSCont.opts.fit.GAP.diff2 = [];
            end
            if ~isfield(MRSCont.opts.fit.GAP, 'diff3')
                MRSCont.opts.fit.GAP.diff3 = [];
            end
        end
        if ~isfield(MRSCont.opts.fit, 'MeanMM')
            MRSCont.opts.fit.MeanMM = 0;
        end
    case 'HERCULES'
        MRSCont.flags.isHERCULES    = 1;
        if exist('editTarget','var')
            MRSCont.opts.editTarget             = editTarget;
        else
            MRSCont.opts.editTarget = {'GABA','GSH'};
        end
        MRSCont.opts.fit.style = opts.fit.style;
        if strcmp(opts.fit.style, 'Concatenated')
            fprintf('Fitting style was changed to Separate, because concatenated modeling is still under development.\n');
            MRSCont.opts.fit.style = 'Separate';
        end
        if isfield(opts.fit, 'coMM3')
            MRSCont.opts.fit.coMM3 = opts.fit.coMM3;
            if isfield(opts.fit, 'FWHMcoMM3')
                MRSCont.opts.fit.FWHMcoMM3 = opts.fit.FWHMcoMM3;
            end
        else
            MRSCont.opts.fit.coMM3 = 'fixedGauss';
            MRSCont.opts.fit.FWHMcoMM3 = 14;
        end
        if ~isfield(MRSCont.opts.fit, 'GAP')
            MRSCont.opts.fit.GAP.sum = [];
            MRSCont.opts.fit.GAP.diff1 = [];
            MRSCont.opts.fit.GAP.diff2 = [];
        else if ~isfield(MRSCont.opts.fit.GAP, 'sum')
                MRSCont.opts.fit.GAP.sum = [];
            end
            if ~isfield(MRSCont.opts.fit.GAP, 'diff1')
                MRSCont.opts.fit.GAP.diff1 = [];
            end
            if ~isfield(MRSCont.opts.fit.GAP, 'diff2')
                MRSCont.opts.fit.GAP.diff2 = [];
            end
        end
        if ~isfield(MRSCont.opts.fit, 'MeanMM')
            MRSCont.opts.fit.MeanMM = 0;
        end
    otherwise
        error('Invalid job file! seqType must be ''unedited'', ''MEGA'', ''HERMES'', or ''HERCULES''.');
end

% AUtomatically check whehther LCModel is used and turn on RAW export
if strcmp(opts.fit.method, 'LCModel')
    MRSCont.opts.saveLCM = 1;
end

% Parse the data scenario entry
if exist('dataScenario','var')
    switch dataScenario
        case 'invivo'
            MRSCont.flags.isPhantom = 0;
        case 'phantom'
            MRSCont.flags.isPhantom = 1;
            % If phantom data are used, override some default fit options
            MRSCont.opts.fit.bLineKnotSpace = 1.0;
            MRSCont.opts.fit.fitMM          = 0;
        case 'Series'
            MRSCont.flags.isSERIES = 1; 
            MRSCont.flags.isPhantom = 0;
        otherwise
            MRSCont.flags.isPhantom = 0;
            warning('Data scenario must be ''invivo'' or ''phantom'' in the job file, and has been set to ''invivo'' (default).');
    end
else
    MRSCont.flags.isPhantom = 0;
    warning('Data scenario must be ''invivo'' or ''phantom'' in the job file, and has been set to ''invivo'' (default).');
end

% Parse the multi voxel entry
if exist('MultiVoxel','var')
    switch MultiVoxel
        case 'PRIAM'
            MRSCont.flags.isPRIAM = 1;
            MRSCont.SENSE = cell(length(priam_offset));
            for kk = 1:length(priam_offset)
                MRSCont.SENSE{kk}.priam_offset = priam_offset{kk};
                MRSCont.SENSE{kk}.priam_direction = priam_direction{kk};
            end

        case 'MRSI'
            MRSCont.flags.isMRSI = 1;
        otherwise
            warning('Multi voxel must be ''PRIAM'' or ''MRSI''in the job file, and has been set to ''single voxel'' (default).');
    end
    if ~isfield(MRSCont.opts, 'MoCo')
        MRSCont.opts.MoCo.target = 'none';
        MRSCont.opts.MoCo.thresh.thresh = 0.8;
        MRSCont.opts.MoCo.thresh.ph_thresh = 0.9;
        MRSCont.opts.MoCo.thresh.last_resort_thresh = 0.6;    
    end
    if ~isfield(MRSCont.opts.MoCo, 'target')
        MRSCont.opts.MoCo.target = 'full';
    end
    if ~isfield(MRSCont.opts.MoCo, 'thresh')
        MRSCont.opts.MoCo.thresh.thresh = 0.8;
        MRSCont.opts.MoCo.thresh.ph_thresh = 0.9;
        MRSCont.opts.MoCo.thresh.last_resort_thresh = 0.6;  
    end
end

% Parse spectral registration entry
if ~isfield(MRSCont.opts,'SpecReg')
    MRSCont.opts.SpecReg = 'RobSpecReg';
    warning('Spectral registration must be ''RobSpecReg'', ''RestSpecReg'' or ''none''  in the job file, and has been set to ''RobSpecReg'' (default).');
end

if exist('file_stat','var')
    if ~isempty(file_stat)
        MRSCont.file_stat = file_stat;
        MRSCont.flags.hasStatfile = 1;
    else
        MRSCont.flags.hasStatfile = 0;
    end
else
    MRSCont.flags.hasStatfile = 0;
end

%Change deface option for images
if isfield(opts, 'img')
   if isfield(opts.img, 'deface')
        MRSCont.opts.img.deface = opts.img.deface;
   end
end

%%% 5. SAVE FILE/FOLDER NAMES INTO MRSCONT %%%
% Make sure that the mandatory fields (metabolite data; output folder) are
% included in the job file.
if exist('files','var')
    MRSCont.files = cellfun(@strtrim,files,'UniformOutput',false);
else
    error('Invalid job file! A job file needs to contain at least metabolite data in the field ''files''.');
end
if exist('files_mm','var')   %re_mm Adding functionality for MM
    MRSCont.files_mm = cellfun(@strtrim,files_mm,'UniformOutput',false);   %re_mm
else
    MRSCont.files_mm = {};
end   %re_mm
if exist('files_mm_ref','var')
    MRSCont.files_mm_ref = cellfun(@strtrim,files_mm_ref,'UniformOutput',false);
else
    MRSCont.files_mm_ref = {};
end
if exist('files_ref','var')
    MRSCont.files_ref = cellfun(@strtrim,files_ref,'UniformOutput',false);
else
    MRSCont.files_ref = {};
end
if exist('files_w','var')
    MRSCont.files_w = cellfun(@strtrim,files_w,'UniformOutput',false);
else
    MRSCont.files_w = {};
end
if exist('files_nii','var') && ~isempty(files_nii)
    MRSCont.files_nii = cellfun(@strtrim,files_nii,'UniformOutput',false);
else
    MRSCont.files_nii = {};
end
if exist('files_seg','var')
    MRSCont.files_seg = cellfun(@strtrim,files_seg,'UniformOutput',false);
else
    MRSCont.files_seg = {};
end
if exist('files_sense','var')
    MRSCont.files_sense = cellfun(@strtrim,files_sense,'UniformOutput',false);
else
    MRSCont.files_sense = {};
end
if exist('outputFolder','var')
    MRSCont.outputFolder = outputFolder;
else
    error('Invalid job file! A job file needs to specify an output folder.');
end

% Check that each array has an identical number of entries
fieldNames = {'files', 'files_ref', 'files_w','files_mm','files_mm_ref', 'files_nii', 'files_seg', 'files_sense'};
ctr = 0;
for kk = 1:length(fieldNames)
    if isfield(MRSCont, fieldNames{kk})
        if ~isempty(MRSCont.(fieldNames{kk}))
            ctr = ctr + 1;
            whichFieldNames{ctr} = fieldNames{kk};
            numDataSets(ctr)     = length(MRSCont.(fieldNames{kk}));
        end
    end
end

% Check whether the number of entries is identical
isUnique = unique(numDataSets);
if length(isUnique) ~= 1
    msg = fprintf('''%s'' has %i entries, but ', whichFieldNames{1}, numDataSets(1));
    for ll = 2:length(whichFieldNames)
        msg2 = fprintf(' ''%s'' has %i entries, ', whichFieldNames{ll}, numDataSets(ll));
        msg = strcat(msg,msg2);
    end

    error('MyComponent:invalidJob', ['Invalid job file! ' msg '\b.']);
end


%%% 6. SET UP DEFAULT OSPREY COLORMAP AND GUI flag %%%
% Default colormap
colormap.Background     = [255/255 254/255 254/255];
colormap.LightAccent    = [110/255 136/255 164/255];
colormap.Foreground     = [11/255 71/255 111/255];
colormap.Accent         = [254/255 186/255 47/255];
MRSCont.colormap        = colormap;
MRSCont.flags.isGUI     = GUI;

%%% 7. SET FLAGS AND VERSION %%%
MRSCont.flags.didJob        = 1;
MRSCont.loadedJob           = jobFile;
MRSCont.ver.Osp             = ['Osprey ' VersionStruct.Version];


%%% 8. CHECK IF OUTPUT STRUCTURE ALREADY EXISTS IN OUTPUT FOLDER %%%
% Determine output folder
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
[~,jobfilename,jobfileext]  = fileparts(jobFile);
outputFile                  = strrep([jobfilename jobfileext], jobfileext, '.mat');
MRSCont.outputFile          = outputFile;
MRSCont.flags.speedUp        = 0;
if ~GUI
    if exist(fullfile(outputFolder, outputFile), 'file') == 2
        switch debug
            case '00'
                fprintf('Your selected output folder ''%s'' already contains an Osprey output structure: \t%s.\n', outputFolder, outputFile);
                fprintf('You are about to load the job: \t%s.\n', jobFile);
                askOverWriteJob = input('Do you want to overwrite the existing job (y/n)? (Warning: This will delete all data in the data container! You can load the existing container by typing n)   [y]   ','s');
            case '11'
                askOverWriteJob = 'y';
            case '01'
                askOverWriteJob = 'n';
        end
        if isempty(askOverWriteJob)
            askOverWriteJob = 'y';
        end
        if askOverWriteJob=='n' || askOverWriteJob=='N'
            switch debug
                case '00'
                     fprintf('You are about to load the job: \t%s.\nIf new files were added they will be attached, otherwise the MRS Container will just be loaded.\nIf you want to change analysis parameters on an exisiting set of files add the location of an existing MRS Container during the call.\n', jobFile);
                     askloadMRSCont = input('Do you want to load the corresponding MRS Container and attach new files (y/n)? [y]   ','s');
                case '01'
                askloadMRSCont = 'y';
            end    
             if isempty(askloadMRSCont)
                askloadMRSCont = 'y';
             end 
             if askloadMRSCont=='n' || askloadMRSCont=='N'
                disp('Aborted! No new job loaded.');
                return;
                else if askloadMRSCont=='y' || askloadMRSCont=='Y'
                        MRSContNew = MRSCont;
                        load(fullfile(outputFolder, outputFile));
                        kk = 1;
                        while (isempty(setdiff(MRSCont.files(kk),MRSContNew.files(kk))))
                            kk = kk + 1;
                            if kk > length(MRSCont.files)
                                kk = kk - 1; 
                                break
                            end
                        end
                        if ~isempty(setdiff(MRSCont.files(kk),MRSContNew.files(kk)))
                            error('The order of the input files must not change. Please append new files at the end of the list');
                        end
                        if length(MRSCont.files) ~= length(MRSContNew.files)
                             MRSCont.files = MRSContNew.files;
                            if isfield(MRSCont,'files_ref')
                                MRSCont.files_ref = MRSContNew.files_ref;
                            end
                            if isfield(MRSCont,'files_w')
                                MRSCont.files_w = MRSContNew.files_w;
                            end
                            if isfield(MRSCont,'files_nii')
                                MRSCont.files_nii = MRSContNew.files_nii;
                            end
                            if isfield(MRSCont,'files_seg')
                                MRSCont.files_seg = MRSContNew.files_seg;
                            end
                            if isfield(MRSCont,'files_sense')
                                MRSCont.files_sense = MRSContNew.files_sense;
                            end
                        end
                        MRSCont.flags.speedUp        = 1;
                    end
             end
        elseif askOverWriteJob=='y' || askOverWriteJob=='Y'            
            delete(fullfile(outputFolder, 'LogFile.txt'));   
            diary(fullfile(outputFolder, 'LogFile.txt'));
            disp('Continue with loading new job, overwriting existing job.');
            fprintf([jobFile '\n']);
            fprintf(['Timestamp %s ' MRSCont.ver.Osp '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));
        end
    else
        diary(fullfile(outputFolder, 'LogFile.txt'));
        fprintf([jobFile '\n']);
        fprintf(['Timestamp %s ' MRSCont.ver.Osp '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));
    end
else
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    if exist(fullfile(outputFolder, outputFile), 'file') == 2
        switch debug
            case '00'
                askOverWriteJob = questdlg({['Your selected output folder already contains an Osprey output structure: \bf' outputFile], ...
                                              ['\rmYou are about to load the job: \bf', strrep(outputFile,'.mat','.m')], ...
                                              '\rmDo you want to overwrite the existing job?',...
                                              '(Warning\rm: This will delete all data in the data container!)',...
                                              'You can load the existing MRS container by clicking No'}, ...
                'Load jobFile', 'Yes','No', opts);
            case '11'
                askOverWriteJob = 'Yes';
            case '01'
                askOverWriteJob = 'No';
        end
        if strcmp(askOverWriteJob, 'No')
            switch debug
                case '00'
                 askloadMRSCont = questdlg( {['You are about to load the job: \bf', strrep(outputFile,'.mat','.m')], ...
                                             '\rmIf new files were added they will be attached, otherwise the MRS Container will just be loaded.', ...
                                             '\rmIf you want to change analysis parameters on an exisiting set of files add a existing MRS Container in the dialog', ...
                                             '\rmDo you want to load the corresponding MRS Container and attach new files?'}, ...
                                              'Load MRS Container', 'Yes','No',opts);
                case '01'
                    askloadMRSCont = 'Yes';
            end 

             if strcmp(askloadMRSCont, 'No')
                disp('Aborted! No new job loaded.');
                return;
                else if strcmp(askloadMRSCont, 'Yes')
                        MRSContNew = MRSCont;
                        load(fullfile(outputFolder, outputFile));
                        kk = 1;
                        while (isempty(setdiff(MRSCont.files(kk),MRSContNew.files(kk))))
                            kk = kk + 1;
                            if kk > length(MRSCont.files)
                                kk = kk - 1; 
                                break
                            end
                        end
                        if ~isempty(setdiff(MRSCont.files(kk),MRSContNew.files(kk)))
                            error('The order of the input files must not change. Please append new files at the end of the list');
                        end
                        if length(MRSCont.files) ~= length(MRSContNew.files)
                             MRSCont.files = MRSContNew.files;
                            if isfield(MRSCont,'files_ref')
                                MRSCont.files_ref = MRSContNew.files_ref;
                            end
                            if isfield(MRSCont,'files_w')
                                MRSCont.files_w = MRSContNew.files_w;
                            end
                            if isfield(MRSCont,'files_mm')
                                MRSCont.files_mm = MRSContNew.files_mm;
                            end
                            if isfield(MRSCont,'files_mm_ref')
                                MRSCont.files_mm_ref = MRSContNew.files_mm_ref;
                            end
                            if isfield(MRSCont,'files_nii')
                                MRSCont.files_nii = MRSContNew.files_nii;
                            end
                            if isfield(MRSCont,'files_seg')
                                MRSCont.files_seg = MRSContNew.files_seg;
                            end
                            if isfield(MRSCont,'files_sense')
                                MRSCont.files_sense = MRSContNew.files_sense;
                            end
                        end
                        MRSCont.flags.speedUp        = 1;
                    end
             end
            elseif strcmp(askOverWriteJob, 'Yes')
                delete(fullfile(outputFolder, 'LogFile.txt'));
                diary(fullfile(outputFolder, 'LogFile.txt'));
                disp('Continue with loading new job, overwriting existing job.');
                fprintf([jobFile '\n']);
                fprintf(['Timestamp %s ' MRSCont.ver.Osp '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));
        end
    else
    diary(fullfile(outputFolder, 'LogFile.txt'));
    fprintf([jobFile '\n']);
    fprintf(['Timestamp %s ' MRSCont.ver.Osp  '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));
    end
end
MRSCont.nDatasets = size(MRSCont.files,2);
%%% 8. SAVE THE OUTPUT STRUCTURE TO THE OUTPUT FOLDER %%%
% Add a save dummy with the GUI flag turned off
saveMRSCont = MRSCont;
MRSCont.flags.isGUI = 0;
MRSCont.flags.isToolChecked = 0;
save(fullfile(outputFolder, outputFile), 'MRSCont');
MRSCont = saveMRSCont;

% Close any remaining open figures
diary off
close all;

end
