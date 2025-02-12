%io_loadspec_twix.m
%Jamie Near, McGill University 2014.
%Edits from Franck Lamberton, 2017.
%
% USAGE:
% out=io_loadspec_twix(filename);
%
% DESCRIPTION:
% Reads in siemens twix raw data (.dat file) using the mapVBVD.m and
% twix_map_obj.m functions from Philipp Ehses (philipp.ehses@tuebingen.mpg.de).
%
% op_loadspec_twix outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
%
% INPUTS:
% filename   = filename of Siemens twix data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.
% out_w      = Input water reference data (only available for some
%               sequences.  This will be empty for others).

function [out,out_w]=io_loadspec_twix(filename);


%read in the data using the new mapVBVD.  This code has been adapted to
%handle both single RAID files and multi-RAID files.  The vast majority of
%Siemens twix data comes as a single RAID file, but I've encoundered a few
%multi-RAID files, particularly when using VD13D.  The way to distinguish
%them here is that a for a single RAID file, mapVBVD will output a struct,
%whereas for a multi-RAID file, mapVBVD will output a cell array of structs.
%This code assumes that the data of interest is in the last element of the
%cell array (possibly a bad assumption under some circumstances):
twix_obj=mapVBVD(filename);
if isstruct(twix_obj)
    disp('single RAID file detected.');
    RaidLength=1;
elseif iscell(twix_obj)
    disp('multi RAID file detected.');
    RaidLength=length(twix_obj);
    %assume that the data of interest is in the last element of the cell.
    twix_obj=twix_obj{RaidLength};
end
dOut.data=twix_obj.image();
version=twix_obj.image.softwareVersion;
sqzSize=twix_obj.image.sqzSize;
sqzDims=twix_obj.image.sqzDims;

%Check if the twix file is from a VE version
if contains(twix_obj.hdr.Dicom.SoftwareVersions, 'E11')
    version='ve';
end
if contains(twix_obj.hdr.Dicom.SoftwareVersions, 'XA')
    index = strfind(twix_obj.hdr.Dicom.SoftwareVersions,'XA');
    version=twix_obj.hdr.Dicom.SoftwareVersions(index:index+3);
end

%find out what sequence, the data were acquired with.  If this is a
%multi-raid file, then the header may contain multiple instances of
%'tSequenceFileName' for different scans (including a pre-scan).
%Therefore, if multi-raid file, we will need to do a bit of extra digging
%to find the correct sequence name.
sequence=twix_obj.hdr.Config.SequenceFileName;

%Try to find out what sequnece this is:
isSpecial=~isempty(strfind(sequence,'rm_special')) ||...  %Is this Ralf Mekle's SPECIAL sequence?
            ~isempty(strfind(sequence,'lr_special')) || ... %This might be Layla Riemann?
            ~isempty(strfind(sequence,'vq_special'));  %or the CIBM SPECIAL sequence?
isjnSpecial=~isempty(strfind(sequence,'jn_svs_special')) ||...  %or Jamie Near's SPECIAL sequence?
            ~isempty(strfind(sequence,'md_Adiab_Special')) ||... %or Masoumeh Dehghani's Adiabatic SPECIAL sequence?
            ~isempty(strfind(sequence,'md_Special')) ||... %or another version of Masoumeh Dehghani's SPECIAL sequence?
            ~isempty(strfind(sequence,'md_Inv_special')); %or Masoumeh Dehghani's Inversion Recovery SPECIAL sequence?
ishdSPECIAL=~isempty(strfind(sequence,'md_dvox_special')); %Is this Masoumeh Dehghani's hadamard-encoded dual-SPECIAL sequence?
isjnMP=~isempty(strfind(sequence,'jn_MEGA_GABA')); %Is this Jamie Near's MEGA-PRESS sequence?
isjnseq=~isempty(strfind(sequence,'jn_')); %Is this another one of Jamie Near's sequences?
isWIP529=~isempty(strfind(sequence,'edit_529'));%Is this WIP 529 (MEGA-PRESS)?
ismodWIP=((~isempty(strfind(sequence,'\svs_edit'))||...
        ~isempty(strfind(sequence,'\wip_svs_edit'))) && ...
        (isempty(strfind(sequence,'edit_859')) && isempty(strfind(sequence,'univ')))); %Modified WIP
isWIP859=~isempty(strfind(sequence,'edit_859'))||...;%Is this WIP 859 (MEGA-PRESS)?
        ~isempty(strfind(sequence,'WIP_859'));%Is this WIP 859 (MEGA-PRESS)? Other naming
isTLFrei=~isempty(strfind(sequence,'md_svs_edit')) ||... %Is Thomas Lange's MEGA-PRESS sequence
         ~isempty(strfind(sequence,'md_svs_slaser_edit')); %Is Thomas Lange's MEGA-s-LASER sequence
isMinn_eja=~isempty(strfind(sequence,'eja_svs_')); %Is this one of Eddie Auerbach's (CMRR, U Minnesota) sequences?
isMinn_dkd=~isempty(strfind(sequence,'svs_slaser_dkd')) ||...   % ... or Dinesh Deelchand's 2016 sLASER sequence?
       ~isempty(strfind(sequence,'svs_slaserVOI_dkd'))||...% ... or Dinesh Deelchand's 2022 'plug and play' sLASER sequence?
       ~isempty(strfind(sequence,'svs_slaserVOI_dkd2'));%Is this Dinesh Deelchand's (CMRR, U Minnesota) sequence?     
isSiemens=(~isempty(strfind(sequence,'svs_se')) ||... %Is this the Siemens PRESS seqeunce?
            ~isempty(strfind(sequence,'svs_st'))) && ... % or the Siemens STEAM sequence?
            isempty(strfind(sequence,'eja_svs'));    %And make sure it's not 'eja_svs_steam'.
isUniversal = ~isempty(strfind(sequence,'univ')) ||...                  % Is JHU universal editing sequence
              ~isempty(strfind(sequence,'smm_svs_herc')) ||...          % Is Pavi's HERCULES sequence
              ~isempty(strfind(sequence,'svs_herc_gls')) ||...          % Is Gize's HERCULES sequence
              ~isempty(strfind(sequence,'svs_se_herc_gls_v2')) || ...   % Is Gize's HERCULES sequence v2
              ~isempty(strfind(sequence,'smm_svs_herm'));               % Is Pavi's HERMES sequence
isDondersMRSfMRI = contains(sequence,'moco_nav_set'); %Is combined fMRI-MRS sequence implemented at Donders Institute NL
isConnectom = contains(twix_obj.hdr.Dicom.ManufacturersModelName,'Connectom'); %Is from Connectom scanner (Apparently svs_se Dims are not as expected for vd)

%Make a pulse sequence identifier for the header (out.seq);
if isSpecial || isjnSpecial
    seq = 'SPECIAL';
elseif isUniversal
    if twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{8} == 1
        seq = 'HERMES';
    else if twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{8} == 3
            seq = 'HERCULES';
        else if isempty(twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{8})
                seq = 'MEGAPRESS';
            else
                seq = 'HERCULES';
            end
        end
    end
elseif isMinn_eja || isMinn_dkd
    if ~isempty(strfind(sequence,'mslaser'))
        seq = 'MEGASLASER';
    else
        if ~isempty(strfind(sequence,'slaser'))
            seq = 'SLASER';
        elseif ~isempty(strfind(sequence,'steam'))
            seq = 'STEAM';
        elseif ~isempty(strfind(sequence,'mpress'))
            seq = 'MEGAPRESS';
        else
            seq = 'PRESS';
        end
    end
elseif isTLFrei && ~isempty(strfind(sequence,'md_svs_slaser_edit'))
    seq = 'MEGASLASER';
elseif isjnMP || isWIP529 || isWIP859 || isTLFrei || ismodWIP
    seq = 'MEGAPRESS';
elseif isSiemens
    if ~isempty(strfind(sequence,'svs_st'))
        seq = 'STEAM';
    elseif ~isempty(strfind(sequence,'svs_se'))
        seq = 'PRESS';
    end
end
if isDondersMRSfMRI
    seq = 'SLASER_D';
end

% GO 10/2022:
% If none of the above match, we've had sequences default to HERMES and a
% bunch of people reported errors... I'll change this to 'PRESS' for now.
if ~exist('seq')
    seq = 'PRESS';
    %seq = 'HERMES';
end

%If this is the SPECIAL sequence, it probably contains both inversion-on
%and inversion-off subspectra on a single dimension, unless it is the VB
%version of Jamie Near's SPECIAL sequence, in which case the subspecs are
%already stored on separate dimensions.
%Both Ralf Mekle's SPECIAL and the VD-VE version of Jamie Near's SPECIAL sequence
%do not store the subspectra along a separate dimension of the data array,
%so we will separate them artifically:
%25 Oct 2018: Due to a recent change, the VE version of Jamie Near's MEGA-PRESS
%sequence also falls into this category.
if isSpecial ||... %Catches Ralf Mekle's and CIBM version of the SPECIAL sequence
        ((strcmp(version,'vd') || strcmp(version,'ve') || contains(version,'XA')) && isjnSpecial) ||... %and the VD/VE versions of Jamie Near's SPECIAL sequence
        ((strcmp(version,'vd') || strcmp(version,'ve') || contains(version,'XA')) && isjnMP);  %and the VD/VE versions of Jamie Near's MEGA-PRESS sequence
    squeezedData=squeeze(dOut.data);
    if twix_obj.image.NCol>1 && twix_obj.image.NCha>1
        data(:,:,:,1)=squeezedData(:,:,[1:2:end-1]);
        data(:,:,:,2)=squeezedData(:,:,[2:2:end]);
        sqzSize=[sqzSize(1) sqzSize(2) sqzSize(3)/2 2];
    elseif twix_obj.NCol>1 && twixObj.image.NCha==1
        data(:,:,1)=squeezedData(:,[1:2:end-1]);
        data(:,:,2)=squeezedData(:,[2:2:end]);
        sqzSize=[sqzSize(1) sqzSize(2)/2 2];
    end
    if isjnseq
        sqzDims{end+1}='Set';
    else
        sqzDims{end+1}='Ida';
    end
elseif ishdSPECIAL %For Masoumeh Dehghani's hadamard-encoded dual-voxel SPECIAL sequence:
    squeezedData=squeeze(dOut.data);
    if twix_obj.image.NCol>1 && twix_obj.image.NCha>1
        data(:,:,:,1)=squeezedData(:,:,[1:4:end-3]);
        data(:,:,:,2)=squeezedData(:,:,[2:4:end-2]);
        data(:,:,:,3)=squeezedData(:,:,[3:4:end-1]);
        data(:,:,:,4)=squeezedData(:,:,[4:4:end]);
        sqzSize=[sqzSize(1) sqzSize(2) sqzSize(3)/4 4];
    elseif twix_obj.image.NCol>1 && twix_obj.image.NCha==1
        data(:,:,1)=squeezedData(:,[1:4:end-3]);
        data(:,:,2)=squeezedData(:,[2:4:end-2]);
        data(:,:,3)=squeezedData(:,[3:4:end-1]);
        data(:,:,4)=squeezedData(:,[4:4:end]);
        sqzSize=[sqzSize(1) sqzSize(2)/4 4];
    end
    if isjnseq
        sqzDims{end+1}='Set';
    else
        sqzDims{end+1}='Ida';
    end
else
    data=dOut.data;
end

%Squeeze the data to remove singleton dims
fids=double(squeeze(data));

%noticed that in the Siemens PRESS and STEAM sequences, there is sometimes
%an extra dimension containing unwanted reference scans or something.  Remove them here.
if isSiemens && (strcmp(version,'vd') || strcmp(version,'ve') || contains(version,'XA')) && strcmp(sqzDims{end},'Phs')
    sqzDims=sqzDims(1:end-1);
    sqzSize=sqzSize(1:end-1);
    if ndims(fids)==4
        fids=fids(:,:,:,2);
        fids=squeeze(fids);
    elseif ndims(fids)==3
        fids=fids(:,:,2);
        fids=squeeze(fids);
    elseif ndims(fids)==2
        fids=fids(:,2);
        fids=squeeze(fids);
    end
end


%Find the magnetic field strength:
Bo=twix_obj.hdr.Dicom.flMagneticFieldStrength;

%Find the number of averages:
Naverages=twix_obj.hdr.Meas.Averages;
if isempty(Naverages)
    Naverages=twix_obj.hdr.MeasYaps.lAverages;
end

%Find out if multiple coil elements were used:
Ncoils=twix_obj.hdr.Meas.iMaxNoOfRxChannels;

%Find the TE:
if ~isMinn_dkd
    TE = twix_obj.hdr.MeasYaps.alTE{1};  %Franck Lamberton
else
    TE = twix_obj.hdr.MeasYaps.alTE{1} + twix_obj.hdr.MeasYaps.alTE{2} + twix_obj.hdr.MeasYaps.alTE{3};
end

%Find the TR:
TR = twix_obj.hdr.MeasYaps.alTR{1};  %Franck Lamberton

wRefs=false;  %Flag to identify if there are automatic water reference scans acquired.  There are none by default.

%Noticed that in Dinesh Deelchand's CMRR sLASER seuqence
%(svs_slaserVOI_dkd2) contains some reference scans, which are acquired at
%the very end.  We will save these and store them as a separate water 
%reference structure (out_w).  
if isMinn_dkd
    %If nRefs is non-zero, then there are automatic water reference scans.
    %These will be saved as "outw".
    nRefs=twix_obj.hdr.MeasYaps.sSpecPara.lAutoRefScanNo;
    if nRefs==1
        nRefs=0;%Becuase Nrefs seems to have a default value of 1 even if there are no reference scans. 
    end
    if nRefs>0
        wRefs=true;
    end
   
    if ndims(fids)==4
        fids_w=fids(:,:,:,end-(nRefs-1):end);
        fids=fids(:,:,:,end-(nRefs+Naverages-1):end-nRefs);
    elseif ndims(fids)==3
        fids_w=fids(:,:,end-(nRefs-1):end);
        fids=fids(:,:,end-(nRefs+Naverages-1):end-nRefs);
    elseif ndims(fids)==2
        fids_w=fids(:,end-(nRefs-1):end);
        fids=fids(:,end-(nRefs+Naverages-1):end-nRefs);
    end
end

% Product Siemens PRESS/STEAM can also include water reference scans. 
if isSiemens && (matches(seq,'PRESS')||matches(seq,'STEAM')) && length(sqzSize)>3
    wRefs=true;
    
    % First, separate out the reference scans:
    Ind = find(matches(sqzDims,'Phs'));
    if ~(sqzSize(Ind)==2)
        error('Expecting 2 entries in the "Phs" dimension, we found %i!',Ind);
    end
    fids_w = squeeze(GetSlice(fids,Ind,1));
    fids   = squeeze(GetSlice(fids,Ind,2));
    sqzDims(Ind) = [];
    sqzSize(Ind) = [];
    
    % Then, concat the averages and repetitions (if they exist)
    if sqzSize(matches(sqzDims,'Rep'))>1
        DimAv = find(matches(sqzDims,'Ave'));
        DimRep = find(matches(sqzDims,'Rep'));
        
        % Reshape the fids:
        NewSize = sqzSize;
        NewSize(DimAv) = NewSize(DimAv) * NewSize(DimRep);
        NewSize(DimRep) = [];
        fids = reshape(fids,NewSize);

        sqzSize = size(fids);
        sqzDims(DimRep)=[];        
        
        % Reshape the water fids, throwing out the zero entries
        NZ_Avs = find(any(~(squeeze(sum(abs(fids_w),[1,2]))==0),DimRep-2)); % Finds water fids in "averages" which are non-zero
        fids_w =fids_w(:,:,NZ_Avs,:);
        Sz = size(fids_w);
        fids_w = reshape(fids_w,Sz(1),Sz(2),[]);
    end

end

% Extract voxel dimensions
if (strcmp(version,'vd') || strcmp(version,'vb') || contains(version,'XA'))
    TwixHeader.VoI_RoFOV     = twix_obj.hdr.Config.VoI_RoFOV; % Voxel size in readout direction [mm]
    TwixHeader.VoI_PeFOV     = twix_obj.hdr.Config.VoI_PeFOV; % Voxel size in phase encoding direction [mm]
    TwixHeader.VoIThickness  = twix_obj.hdr.Config.VoI_SliceThickness; % Voxel size in slice selection direction [mm]
    TwixHeader.PosCor         = twix_obj.hdr.Config.VoI_Position_Cor; % Coronal coordinate of voxel [mm]
    TwixHeader.PosSag         = twix_obj.hdr.Config.VoI_Position_Sag; % Sagittal coordinate of voxel [mm]
    TwixHeader.PosTra         = twix_obj.hdr.Config.VoI_Position_Tra; % Transversal coordinate of voxel [mm]
    TwixHeader.VoI_InPlaneRot = twix_obj.hdr.Config.VoI_InPlaneRotAngle; % Voxel rotation in plane
    TwixHeader.NormCor        = twix_obj.hdr.Config.VoI_Normal_Cor; % Coronal component of normal vector of voxel
    TwixHeader.NormSag        = twix_obj.hdr.Config.VoI_Normal_Sag; % Sagittal component of normal vector of voxel
    TwixHeader.NormTra        = twix_obj.hdr.Config.VoI_Normal_Tra; % Transversal component of normal vector of voxel
else
    TwixHeader.VoI_RoFOV     = twix_obj.hdr.Spice.VoiReadoutFOV; % Voxel size in readout direction [mm]
    TwixHeader.VoI_PeFOV     = twix_obj.hdr.Spice.VoiPhaseFOV; % Voxel size in phase encoding direction [mm]
    TwixHeader.VoIThickness  = twix_obj.hdr.Spice.VoiThickness; % Voxel size in slice selection direction [mm]
    TwixHeader.PosCor         = twix_obj.hdr.Spice.VoiPositionCor; % Coronal coordinate of voxel [mm]
    TwixHeader.PosSag         = twix_obj.hdr.Spice.VoiPositionSag; % Sagittal coordinate of voxel [mm]
    TwixHeader.PosTra         = twix_obj.hdr.Spice.VoiPositionTra; % Transversal coordinate of voxel [mm]
    TwixHeader.VoI_InPlaneRot = twix_obj.hdr.Spice.VoiInPlaneRot; % Voxel rotation in plane
    TwixHeader.NormCor        = twix_obj.hdr.Spice.VoiNormalCor; % Coronal component of normal vector of voxel
    TwixHeader.NormSag        = twix_obj.hdr.Spice.VoiNormalSag; % Sagittal component of normal vector of voxel
    TwixHeader.NormTra        = twix_obj.hdr.Spice.VoiNormalTra; % Transversal component of normal vector of voxel
end
TwixHeader.TablePosSag    = twix_obj.hdr.Dicom.lGlobalTablePosSag; % Sagittal table position [mm]
TwixHeader.TablePosCor    = twix_obj.hdr.Dicom.lGlobalTablePosCor; % Coronal table position [mm]
TwixHeader.TablePosTra    = twix_obj.hdr.Dicom.lGlobalTablePosTra; % Transversal table position [mm]
% If a parameter is set to zero (e.g. if no voxel rotation is
% performed), the respective field is left empty in the TWIX file. This
% case needs to be intercepted. Setting to the minimum possible value.
VoI_Params = {'VoI_InPlaneRot','VoI_RoFOV','VoI_PeFOV','VoIThickness','NormCor','NormSag','NormTra', ...
              'PosCor','PosSag','PosTra','TablePosSag','TablePosCor','TablePosTra'};
for pp = 1:length(VoI_Params)
    if isempty(TwixHeader.(VoI_Params{pp}))
        TwixHeader.(VoI_Params{pp}) = realmin('double');
    end
end
geometry.size.VoI_RoFOV     = TwixHeader.VoI_RoFOV; % Voxel size in readout direction [mm]
geometry.size.VoI_PeFOV     = TwixHeader.VoI_PeFOV; % Voxel size in phase encoding direction [mm]
geometry.size.VoIThickness  = TwixHeader.VoIThickness; % Voxel size in slice selection direction [mm]
geometry.pos.PosCor         = TwixHeader.PosCor; % Coronal coordinate of voxel [mm]
geometry.pos.PosSag         = TwixHeader.PosSag; % Sagittal coordinate of voxel [mm]
geometry.pos.PosTra         = TwixHeader.PosTra; % Transversal coordinate of voxel [mm]
geometry.pos.TablePosSag    = TwixHeader.TablePosSag; % Sagittal table position [mm]
geometry.pos.TablePosCor    = TwixHeader.TablePosCor; % Coronal table position [mm]
geometry.pos.TablePosTra    = TwixHeader.TablePosTra; % Transversal table position [mm]
geometry.rot.VoI_InPlaneRot = TwixHeader.VoI_InPlaneRot; % Voxel rotation in plane
geometry.rot.NormCor        = TwixHeader.NormCor; % Coronal component of normal vector of voxel
geometry.rot.NormSag        = TwixHeader.NormSag; % Sagittal component of normal vector of voxel
geometry.rot.NormTra        = TwixHeader.NormTra; % Transversal component of normal vector of voxel

%Now begin indexing the dimensions of the data array. ie. create the dims
%structure, which specifies which dimensions of the data array are being
%used to hold the time-domain data, the multiple coil channels, the
%average, the sub-spectra, and any additional dimensions.
sqzDims_update=sqzDims;
dimsToIndex=[1:length(sqzDims)];

%First index the dimension of the time-domain data
dims.t=find(strcmp(sqzDims,'Col'));
if ~isempty(dims.t)
    %remove the time dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.t);
else
    dims.t=0;
    error('ERROR:  Spectrom contains no time domain information!!');
end

%Now index the dimension of the coil channels
dims.coils=find(strcmp(sqzDims,'Cha'));
if ~isempty(dims.coils)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.coils);
else
    dims.coils=0;
end

%Now index the dimension of the averages
if strcmp(version,'vd') || strcmp(version,'ve') || contains(version,'XA')
    if isMinn_eja || isMinn_dkd || isConnectom
        dims.averages=find(strcmp(sqzDims,'Set'));
    else
        % GO 6/2024: Encountered a VD13 Siemens product sequence that had
        % transients in Set, not in Ave... 
        dims.averages=find(strcmp(sqzDims,'Ave'));
        if isempty(dims.averages)
            dims.averages=find(strcmp(sqzDims,'Set'));
        end
    end
else
    dims.averages=find(strcmp(sqzDims,'Set'));
end
if ~isempty(dims.averages)
    %remove the averages dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
else
    %If no Averages dimension was found, then check for a "Repetitions"
    %dimension.  If that is found, store it under "averages".  If both
    %"Averages" and "Repetitions" dimensions are found, "Repetitions" will
    %be indexed under "Extras", since "Repetitions is not currently an
    %option in FID-A.
    dims.averages=find(strcmp(sqzDims,'Rep'));
    if ~isempty(dims.averages)
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
    else
        %If neither an "Averages" or a "Repetitions" dimension is found,
        %then set the FID-A "Averages" dimension to zero.
        dims.averages=0;
    end
end

%Now we have indexed the dimensions containing the timepoints, the coil
%channels, and the averages.  As we indexed each dimension, we removed the
%corresponding index from the dimsToIndex vector.  At this point, if there
%are any values left in the dimsToIndex vector, then there must be some
%additional dimensions that need indexing.  We assume that if sub-spectra exist,
%then these must be indexed in either the 'Ida' dimension (for all Jamie
%Near's VB-version pulse sequences), the 'Set' dimension (for all Jamie
%Near's VD/VE-version pulse sequences), the 'Eco' dimension (for the WIP
%529 MEGA-PRESS sequence or the Minnesota MEGA-PRESS sequence), or the 'Ide'
% dimension (for the WIP 859 MEGA-PRESS sequence).
if ~isempty(dimsToIndex)
    %Now index the dimension of the sub-spectra
    if isjnseq  || isSpecial
        if strcmp(version,'vd') || strcmp(version,'ve')
            dims.subSpecs=find(strcmp(sqzDims,'Set'));
        else
            dims.subSpecs=find(strcmp(sqzDims,'Ida'));
        end
    elseif isWIP529 || isMinn_eja || isMinn_dkd
        dims.subSpecs=find(strcmp(sqzDims,'Eco'));
    elseif ismodWIP
        if strcmp(version,'vd') || strcmp(version,'ve') || strcmp(version,'vb')
            dims.subSpecs=find(strcmp(sqzDims,'Eco'));
        else
            dims.subSpecs=find(strcmp(sqzDims,'Ide'));
        end
    elseif isWIP859
        dims.subSpecs=find(strcmp(sqzDims,'Ide'));
    else
        dims.subSpecs=dimsToIndex(1);
    end
    if ~isempty(dims.subSpecs)
        %remove the sub-spectra dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.subSpecs);
    else
        dims.subSpecs=0;
    end
else
    dims.subSpecs=0;
end

%And if any further dimensions exist after indexing the sub-spectra, call
%these the 'extras' dimension.
if ~isempty(dimsToIndex)
    %Now index the 'extras' dimension
    dims.extras=dimsToIndex(1);
    if ~isempty(dims.extras)
        %remove the extras dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.extras);
    else
        dims.extras=0;
    end
else
    dims.extras=0;
end

%Now that we've indexed the dimensions of the data array, we now need to
%permute it so that the order of the dimensions is standardized:  we want
%the order to be as follows:
%   1) time domain data.
%   2) coils.
%   3) averages.
%   4) subSpecs.
%   5) extras.
if length(sqzDims)==5
    fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs dims.extras]);
    dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.extras=5;
    if wRefs
        fids_w=permute(fids_w,[dims.t dims.coils dims.averages dims.subSpecs dims.extras]);
    end
elseif length(sqzDims)==4
    if dims.extras==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.extras=0;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.coils dims.averages dims.subSpecs]);
        end
    elseif dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.extras=4;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.coils dims.averages dims.extras]);
        end
    elseif dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.extras=4;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.coils dims.subSpecs dims.extras]);
        end
    elseif dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.extras=4;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.averages dims.subSpecs dims.extras]);
        end
    end
elseif length(sqzDims)==3
    if dims.extras==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.extras=0;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.coils dims.averages]);
        end
    elseif dims.extras==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.extras=0;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.coils dims.subSpecs]);
        end
    elseif dims.extras==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.extras=0;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.averages dims.subSpecs]);
        end
    end
elseif length(sqzDims)==2
    if dims.extras==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.extras=0;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.coils]);
        end
    elseif dims.extras==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.extras=0;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.averages]);
        end
    elseif dims.extras==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.extras=0;
        if wRefs
            fids_w=permute(fids_w,[dims.t dims.subSpecs]);
        end
    end
elseif length(sqzDims)==1
    fids=permute(fids,[dims.t]);
    dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.extras=0;
    if wRefs
        fids_w=permute(fids_w,[dims.t]);
    end
end

%Now reorder the fids for the Universal MEGA implementation
if strcmp(seq,'MEGAPRESS') && isUniversal
    % Can only do this if we have more than one row, which might not be the
    % case (e.g. for water)
    if size(fids,3) > 1
        fids_A = fids(:,:,1:2:end);
        fids_B = fids(:,:,2:2:end);
        fids = cat(4,fids_A,fids_B);
        dims.subSpecs=4;
    else
        dims.subSpecs=0;
    end
end

%Now reorder the fids for the Universal HERMES/HERCULES implementation
if strcmp(seq,'HERMES') || strcmp(seq,'HERCULES')
    % Can only do this if we have more than one row, which might not be the
    % case (e.g. for water)
    if size(fids,3) > 1
        fids_A = fids(:,:,1:4:end);
        fids_B = fids(:,:,2:4:end);
        fids_C = fids(:,:,3:4:end);
        fids_D = fids(:,:,4:4:end);
        fids = cat(4,fids_A,fids_B,fids_C,fids_D);
        if size(fids,3) > 1
            dims.subSpecs=4;
        else
            fids = squeeze(fids);
            dims.subSpecs=3;
            dims.averages=0;
        end
    else
        dims.subSpecs=0;
    end
end

%Now get the size of the data array:
sz=size(fids);
if wRefs
    sz_w=size(fids_w);
end

%Now take fft of time domain to get fid:
specs=fftshift(fft(fids,[],dims.t),dims.t);
if wRefs
    specs_w=fftshift(ifft(fids_w,[],dims.t),dims.t);
end

%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;  %Franck Lamberton
spectralwidth=1/dwelltime;

%Get TxFrq
if strcmp(version,'ve')
    txfrq=twix_obj.hdr.Meas.lFrequency  ;
else
    try
        txfrq=twix_obj.hdr.Meas.Frequency;
    catch
        txfrq=twix_obj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;
    end
    if isempty(txfrq)
        txfrq=twix_obj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;
    end
end

%Get Date
%date = getfield(regexp(twix_obj.hdr.MeasYaps.tReferenceImage0, ...
%'^".*\.(?<DATE>\d{8})\d*"$', 'names'), 'DATE');  %Franck Lamberton

date=''; %The above code for extracting the date from the header
         %was causing problems.  Since date is not critical
         %for almost any applications, removing it now to be fixed at a
         %later date.

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired
%averages in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    if dims.averages~=0
        averages=sz(dims.averages)*sz(dims.subSpecs);
        rawAverages=averages;
        if wRefs
            averages_w=sz_w(dims.averages)*sz_w(dims.subSpecs);
            rawAverages_w=averages_w;
        end
    else
        averages=sz(dims.subSpecs);
        rawAverages=1;
        if wRefs
            averages_w=sz_w(dims.subSpecs);
            rawAverages_w=1;
        end
    end
else
    if dims.averages~=0
        averages=sz(dims.averages);
        rawAverages=averages;
        if wRefs
            averages_w=sz_w(dims.averages);
            rawAverages_w=averages_w;
        end
    else
        averages=1;
        rawAverages=1;
        if wRefs
            averages_w=1;
            rawAverages_w=1;
        end
    end
end

%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired
%subspectra in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end

%Find the number of points acquired before the echo so that this
%information can be stored in the .pointsToLeftshfit field of the data
%structure.  Depending on the pulse sequence used to acquire the data, the
%header location of this parameter is different.  For product PRESS
%seqeunces, the value is located in twix_obj.image.freeParam(1).  For WIP
%sequences, the value is located in twix_obj.image.cutOff(1,1).  For CMRR
%sequences, the value is located in twix_obj.image.iceParam(5,1).  Special
%thanks to Georg Oeltzschner for decoding all of this and sharing the
%information with me:

if isWIP529 || isWIP859
    leftshift = twix_obj.image.cutOff(1,1);
elseif isSiemens && (~(isMinn_eja || isMinn_dkd) && ~isConnectom)
    if ~strcmp(version,'ve') && ~contains(version,'XA')
        leftshift = twix_obj.image.freeParam(1);
    else
       leftshift = twix_obj.image.iceParam(5,1);
    end
elseif isMinn_eja || isMinn_dkd || isConnectom || isDondersMRSfMRI
    try
        leftshift = twix_obj.image.iceParam(5,1);
    catch
        leftshift = twix_obj.image.freeParam(1);
    end
elseif isTLFrei && strcmp(seq,'MEGASLASER')
    leftshift = twix_obj.image.iceParam(5,1);
else
    leftshift = twix_obj.image.freeParam(1);
end

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=f/(Bo*42.577);

% Siemens data assumes the center frequency to be 4.7 ppm:
centerFreq = 4.7;
ppm=ppm + centerFreq;

t=[0:dwelltime:(sz(1)-1)*dwelltime];


%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;
out.t=t;
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=seq;
out.te=TE/1000;
out.tr=TR/1000;
out.pointsToLeftshift=leftshift;
out.centerFreq = centerFreq;
out.geometry = geometry;
if isfield(twix_obj.hdr.Config,'Nucleus')
    out.nucleus = twix_obj.hdr.Config.Nucleus  ;
end
if isfield(twix_obj.hdr.Dicom,'SoftwareVersions')
    out.software = [version ' ' twix_obj.hdr.Dicom.SoftwareVersions];
else
    out.software = version;
end

%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end
% Add info for niiwrite
out.PatientPosition = twix_obj.hdr.Config.PatientPosition;
out.Manufacturer = 'Siemens';
[~,filename,ext] = fileparts(filename);
out.OriginalFile = [filename ext];
% Sequence flags
out.flags.isUnEdited = 0;
out.flags.isMEGA = 0;
out.flags.isHERMES = 0;
out.flags.isHERCULES = 0;
out.flags.isPRIAM = 0;
out.flags.isMRSI = 0;
if strcmp(seq,'PRESS') || strcmp(seq,'STEAM') || strcmp(seq,'SLASER')
    out.flags.isUnEdited = 1;
end
if contains(seq,'MEGA')
    out.flags.isMEGA = 1;
end
if strcmp(seq,'HERMES')
    out.flags.isHERMES = 1;
end
if strcmp(seq,'HERCULES')
    out.flags.isHERCULES = 1;
end

if wRefs
    %FILLING IN DATA STRUCTURE
    out_w.fids=fids_w;
    out_w.specs=specs_w;
    out_w.sz=sz_w;
    out_w.ppm=ppm;
    out_w.t=t;
    out_w.spectralwidth=spectralwidth;
    out_w.dwelltime=dwelltime;
    out_w.txfrq=txfrq;
    out_w.date=date;
    out_w.dims=dims;
    out_w.Bo=Bo;
    out_w.averages=averages_w;
    out_w.rawAverages=rawAverages_w;
    out_w.subspecs=subspecs;
    out_w.rawSubspecs=rawSubspecs;
    out_w.seq=seq;
    out_w.te=TE/1000;
    out_w.tr=TR/1000;
    out_w.pointsToLeftshift=leftshift;
    out_w.centerFreq = centerFreq;
    out_w.geometry = geometry;
    if isfield(twix_obj.hdr.Config,'Nucleus')
        out_w.nucleus = twix_obj.hdr.Config.Nucleus  ;
    end
    if isfield(twix_obj.hdr.Dicom,'SoftwareVersions')
        out_w.software = [version ' ' twix_obj.hdr.Dicom.SoftwareVersions];
    else
        out_w.software = version;
    end


    %FILLING IN THE FLAGS
    out_w.flags.writtentostruct=1;
    out_w.flags.gotparams=1;
    out_w.flags.leftshifted=0;
    out_w.flags.filtered=0;
    out_w.flags.zeropadded=0;
    out_w.flags.freqcorrected=0;
    out_w.flags.phasecorrected=0;
    out_w.flags.averaged=0;
    out_w.flags.addedrcvrs=0;
    out_w.flags.subtracted=0;
    out_w.flags.writtentotext=0;
    out_w.flags.downsampled=0;
    if out_w.dims.subSpecs==0
        out_w.flags.isFourSteps=0;
    else
        out_w.flags.isFourSteps=(out_w.sz(out_wop_pl.dims.subSpecs)==4);
    end
    % Add info for niiwrite
    out_w.PatientPosition = twix_obj.hdr.Config.PatientPosition;
    out_w.Manufacturer = 'Siemens';
    [~,filename,ext] = fileparts(filename);
    out_w.OriginalFile = [filename ext];
    % Sequence flags
    out_w.flags.isUnEdited = 0;
    out_w.flags.isMEGA = 0;
    out_w.flags.isHERMES = 0;
    out_w.flags.isHERCULES = 0;
    out_w.flags.isPRIAM = 0;
    out_w.flags.isMRSI = 0;
    if strcmp(seq,'PRESS') || strcmp(seq,'STEAM') || strcmp(seq,'SLASER')
        out_w.flags.isUnEdited = 1;
    end
    if contains(seq,'MEGA')
        out_w.flags.isMEGA = 1;
    end
    if strcmp(seq,'HERMES')
        out_w.flags.isHERMES = 1;
    end
    if strcmp(seq,'HERCULES')
        out_w.flags.isHERCULES = 1;
    end
else
    %No integrated water reference data found.  Returning empty struct for out_w:
    disp('No integrated water reference data found.  Returning empty field');
    out_w=[];
end

end
%DONE

function[SlicedArray] = GetSlice(Array, Dim, Idx)
%% function[SlicedArray] = GetSlice(Array, Dim, Idx)
%
% Description: Takes a slice of multi-dimensional array in dimension "Dim"
% at the specified index: Idx.
%
% Input:     Array = Aray to be sliced 
%            Dim = Dimension in which to slice
%            Idx = Index (or indices) to extract from dimension "Dim"
% Output:    SlicedArray = Subset of array
%
% Example usage:
%   For a 7D array, "Data":
%
%   SlicedData = GetSlice(Data, 5, 1:5);
%       returns 7D array but with only the 1st 5 entries in 5th dim
%
% C.W. Davies-Jenkins, Johns Hopkins University 2024

Specification = repmat({':'}, 1, ndims(Array)); % Create an vector whose length matches the dimensionality of the Array
Specification{Dim} = Idx; % At the specified Dim, add the indices we want to extract

SlicedArray = Array(Specification{:}); % Apply the slicing!

end
