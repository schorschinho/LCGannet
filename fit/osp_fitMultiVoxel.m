function [outMRSCont] = osp_fitMultiVoxel(MRSCont)    
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
metFitTime = tic;
outMRSCont= MRSCont;
fitMRSCont = MRSCont;
reverseStr = '';
%% Get infos to set up a loop to process all voxels
if MRSCont.flags.isPRIAM == 1
    XVox = MRSCont.raw{1}.nXvoxels;
else if MRSCont.flags.isMRSI == 1
        XVox = MRSCont.raw{1}.nXvoxels;
        YVox = MRSCont.raw{1}.nYvoxels;
        ZVox = MRSCont.raw{1}.nZvoxels;
    end
end
SubSpecNames = fieldnames(fitMRSCont.processed);
NoSubSpec = length(fieldnames(fitMRSCont.processed)); 

if MRSCont.flags.isPRIAM == 1
    for x = 1 : XVox
       for ss = 1 : NoSubSpec % Loop over Subspec 
            for kk = 1 :MRSCont.nDatasets
                    fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},x);                
           end
       end
        for kk = 1 :MRSCont.nDatasets % Loop over scale values
                    fitMRSCont.fit.scale{kk} =  MRSCont.fit.scale{kk}(x);                
        end
        if MRSCont.flags.isUnEdited
            [fitMRSCont] = osp_fitUnEdited(fitMRSCont);
        elseif MRSCont.flags.isMEGA
            [fitMRSCont] = osp_fitMEGA(fitMRSCont);
        elseif MRSCont.flags.isHERMES
            [fitMRSCont] = osp_fitHERMES(fitMRSCont);
        elseif MRSCont.flags.isHERCULES
            % For now, fit HERCULES like HERMES data
            [fitMRSCont] = osp_fitHERCULES(fitMRSCont);
        else
            msg = 'No flag set for sequence type!';
            fprintf(fileID,msg);
            error(msg);
        end

        if x == 1
            outMRSCont.fit = fitMRSCont.fit;
        else
            fields = {'resBasisSet','results'};
            for f = 1 : length(fields)
                if isfield(outMRSCont.fit,fields{f})
                    if iscell(outMRSCont.fit.(fields{f}))
                        %PRIAM data
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
                    else
                        temp = outMRSCont.fit.(fields{f});
                        outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
                        outMRSCont.fit.(fields{f}){1} = temp;
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
                    end            
                end
            end
        end
    end
    time = toc(metFitTime);
    outMRSCont.runtime.FitMet = time;    
elseif MRSCont.flags.isMRSI == 1
    %Fit center of MRSI first for inital guess parameters
   cx = round(XVox/2);
   cy = round(YVox/2);
   cz = round(ZVox/2);

   for ss = 1 : NoSubSpec % Loop over Subspec 
        for kk = 1 :MRSCont.nDatasets
                fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},x);                
       end
        for kk = 1 :MRSCont.nDatasets % Loop over scale values
                    fitMRSCont.fit.scale{kk} =  MRSCont.fit.scale{kk}(x);                
        end
        if MRSCont.flags.isUnEdited
            [fitMRSCont] = osp_fitUnEdited(fitMRSCont);
        elseif MRSCont.flags.isMEGA
            [fitMRSCont] = osp_fitMEGA(fitMRSCont);
        elseif MRSCont.flags.isHERMES
            [fitMRSCont] = osp_fitHERMES(fitMRSCont);
        elseif MRSCont.flags.isHERCULES
            % For now, fit HERCULES like HERMES data
            [fitMRSCont] = osp_fitHERCULES(fitMRSCont);
        else
            msg = 'No flag set for sequence type!';
            fprintf(fileID,msg);
            error(msg);
        end

        if x == 1
            outMRSCont.fit = fitMRSCont.fit;
        else
            fields = {'resBasisSet','results'};
            for f = 1 : length(fields)
                if isfield(outMRSCont.fit,fields{f})
                    if iscell(outMRSCont.fit.(fields{f}))
                        %PRIAM data
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
                    else
                        temp = outMRSCont.fit.(fields{f});
                        outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
                        outMRSCont.fit.(fields{f}){1} = temp;
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
                    end            
                end
            end
        end
    end
elseif MRSCont.flags.isMRSI == 1
    %Fit center of MRSI first for inital guess parameters
   cx = round(XVox/2);
   cy = round(YVox/2);
   cz = round(ZVox/2);
   
   for ss = 1 : NoSubSpec % Loop over Subspec 
       for kk = 1 :MRSCont.nDatasets
           if ZVox <=1
               fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[cx,cy]);  
           else
               fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[cx,cy,cz]); 
           end
       end
   end                       
    msg = sprintf('\nFitting center voxel (%d, %d, %d) of the MRSI dataset...\n', cx, cy, cz);
    fprintf([reverseStr, msg]);

    if MRSCont.flags.isUnEdited
        [fitMRSCont] = osp_fitUnEdited(fitMRSCont);
    elseif MRSCont.flags.isMEGA
        [fitMRSCont] = osp_fitMEGA(fitMRSCont);
    elseif MRSCont.flags.isHERMES
        [fitMRSCont] = osp_fitHERMES(fitMRSCont);
    elseif MRSCont.flags.isHERCULES
        % For now, fit HERCULES like HERMES data
        [fitMRSCont] = osp_fitHERCULES(fitMRSCont);
    else
        msg = 'No flag set for sequence type!';
        fprintf(fileID,msg);
        error(msg);
    end
    outMRSCont.fit = fitMRSCont.fit;
    fields = {'resBasisSet','results'};
    for f = 1 : length(fields)
        if ZVox <=1
            temp = outMRSCont.fit.(fields{f});
            outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
            outMRSCont.fit.(fields{f}){cx,cy} = fitMRSCont.fit.(fields{f});
        else  % 3D MRSI data
            temp = outMRSCont.fit.(fields{f});
            outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
            outMRSCont.fit.(fields{f}){cx,cy,cz} = fitMRSCont.fit.(fields{f});
        end          
    end

   MRSCont.fit.MRSIfitPriors = fitMRSCont.fit;
    % Fit all remaining voxels
    for z = 1 : ZVox 
        for nVox = 2 : 64
           [x,y] = osp_spiral(nVox);
           x = x+cx;
           y = y+cy;
           for ss = 1 : NoSubSpec % Loop over Subspec 
               for kk = 1 :MRSCont.nDatasets
                   if ZVox <=1
                       fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[x,y]);  
                   else
                       fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[x,y,z]); 
                   end
               end
           end
            fitMRSCont.fit =  MRSCont.fit.MRSIfitPriors;  %Load prior results into the struct                        
            msg = sprintf('\nFitting metabolite spectra from kx %d out of %d total x phase steps...\n', x, XVox);
            fprintf([reverseStr, msg]);
            msg = sprintf('\nFitting metabolite spectra from ky  %d out of %d total y phase steps...\n', y, YVox);
            fprintf([reverseStr, msg]);
            msg = sprintf('\nFitting metabolite spectra from slice %d out of %d total slices...\n', z, ZVox);
            fprintf([reverseStr, msg]);
            msg = sprintf('\nFitting metabolite spectra from voxel %d out of %d total voxels...\n', x*y*z, XVox*YVox*ZVox);
            fprintf([reverseStr, msg]);
            if ~((x == cx) && (y == cy) && (z == cz)) % Do not re-analyze the center voxel
                if MRSCont.flags.isUnEdited
                    [fitMRSCont] = osp_fitUnEdited(fitMRSCont);
                elseif MRSCont.flags.isMEGA
                    [fitMRSCont] = osp_fitMEGA(fitMRSCont);
                elseif MRSCont.flags.isHERMES
                    [fitMRSCont] = osp_fitHERMES(fitMRSCont);
                elseif MRSCont.flags.isHERCULES
                    % For now, fit HERCULES like HERMES data
                    [fitMRSCont] = osp_fitHERCULES(fitMRSCont);
                else
                    msg = 'No flag set for sequence type!';
                    fprintf(fileID,msg);
                    error(msg);
                end
                MRSCont.fit.MRSIfitPriors = fitMRSCont.fit; %Store new priors
                fields = {'resBasisSet','results'};
                for f = 1 : length(fields)
                     % 2D MRSI data
                    if ZVox <=1
                        outMRSCont.fit.(fields{f}){x,y} = fitMRSCont.fit.(fields{f});
                    else  % 3D MRSI data
                        outMRSCont.fit.(fields{f}){x,y,z} = fitMRSCont.fit.(fields{f});
                    end          
                end
            end
        end
    end
    time = toc(metFitTime);
    outMRSCont.runtime.FitMet = time;    
    outMRSCont.fit.basisSet = MRSCont.fit.basisSet;
    outMRSCont.fit.scale = MRSCont.fit.scale;
end
end