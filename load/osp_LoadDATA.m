function [MRSCont] = osp_LoadDATA(MRSCont)
%% [MRSCont] = osp_LoadSDAT(MRSCont)
%   Reads Philips DATA/LIST files it is adapted from Gannet.
%
%   Author:
%       Dr.Helge Zoellner (Johns Hopkins University, 2020-10-02)
%       hzoelln2@jhmi.edu
%   
%   Credits:
%       This code uses the function
%       loadRawKspace.m
%       from the excellent "Matlab raw kspace tools" toolbox
%       (Wouter Potters, Academic Medical Center, Amsterdam, NL)
%       https://bitbucket.org/wpotters/matlab-raw-kspace-tools
%
%   History:
%       2020-10-02: First version.


% Close any remaining open figures
close all;
warning('off','all');
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
if MRSCont.flags.hasMM %re_mm adding functionality to load MM data
    if ((length(MRSCont.files_mm) == 1) && (MRSCont.nDatasets>1))   %re_mm seems like specificy one MM file for a batch is also an option to plan to accomodate
        for kk=2:MRSCont.nDatasets %re_mm 
            MRSCont.files_mm{kk} = MRSCont.files_mm{1}; % re_mm allowable to specify one MM file for the whole batch
        end %re_mm 
    end   %re_mm 
    if ((length(MRSCont.files_mm) ~= MRSCont.nDatasets) )   %re_mm 
        msg = 'Number of specified MM files does not match number of specified metabolite files.'; %re_mm 
        fprintf(fileID,msg);
        error(msg);
    end   %re_mm 
end   %re_mm 
if MRSCont.flags.hasRef
    if length(MRSCont.files_ref) ~= MRSCont.nDatasets
        msg = 'Number of specified reference files does not match number of specified metabolite files.'; %re_mm 
        fprintf(fileID,msg);
        error(msg);
    end
end
if MRSCont.flags.hasWater
    if length(MRSCont.files_w) ~= MRSCont.nDatasets
        msg = 'Number of specified water files does not match number of specified metabolite files.'; %re_mm 
        fprintf(fileID,msg);
        error(msg);
    end
end

%% Get the data (loop over all datasets)
refLoadTime = tic;
reverseStr = '';
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
end
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
for kk = 1:MRSCont.nDatasets
msg = sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    fprintf(fileID,[reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    if MRSCont.flags.isGUI        
        set(progressText,'String' ,sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
    end
    
    if ((MRSCont.flags.didLoadData == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~isfield(MRSCont.ver, 'Load') || ~strcmp(MRSCont.ver.Load,MRSCont.ver.CheckLoad))
        
        % Read in the raw metabolite data. Since the Philips DATA loader needs
        % to know the number of sub-spectra (e.g. from spectral editing), the
        % type of sequence needs to be differentiated here already.
        if MRSCont.flags.hasStatfile
            statFile = MRSCont.file_stat;
        else
            statFile = [];
        end
           
        if MRSCont.flags.isUnEdited
            [raw,raw_ref]=io_loadspec_data(MRSCont.files{kk},1,kk,statFile);
        elseif MRSCont.flags.isMEGA
            [raw,raw_ref]=io_loadspec_data(MRSCont.files{kk},2,kk,statFile);
        elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
            [raw,raw_ref]=io_loadspec_data(MRSCont.files{kk},4,kk,statFile);
        end
        MRSCont.raw_uncomb{kk}      = raw;
        
        if ~MRSCont.flags.hasRef && ~isempty(raw_ref)
            MRSCont.raw_ref_uncomb{kk}  = raw_ref;
            MRSCont.flags.hasRef = 1;
        else if MRSCont.flags.hasRef
                [~,raw_ref]=io_loadspec_data(MRSCont.files_ref{kk},1,kk,statFile);
                MRSCont.raw_ref_uncomb{kk}  = raw_ref;
            end
        end
        
        if MRSCont.flags.hasWater
            [~,raw_w]=io_loadspec_data(MRSCont.files_w{kk},1,kk,statFile);
            MRSCont.raw_w_uncomb{kk}    = raw_w;
        end
    end
    
end

fprintf('... done.\n');
time = toc(refLoadTime);
if MRSCont.flags.isGUI        
    set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
    pause(1);
end
fprintf(fileID,'... done.\n Elapsed time %f seconds\n',time);
fclose(fileID);
% Set flag
MRSCont.flags.coilsCombined     = 0;
MRSCont.runtime.Load = time;
end
