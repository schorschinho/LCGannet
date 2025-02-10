%io_loadspec_dicom.m
%Jamie Near, McGill University 2014.
%Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% out=io_loadspec_dicom(source);
%
% DESCRIPTION:
% Loads a DICOM (.dcm, .ima) file into matlab structure format.
%
% INPUTS:
% source     = Name of a folder (or file) with DICOM (.dcm, .ima) data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out=io_loadspec_dicom(source);

% ARC170724 : Source may be a folder, or an individual file. We later determine
% how this should be interpreted.
if ~isfolder(source) && isfile(source)
    [folder,filename,ext] =fileparts(source);
    requested_file=[filename ext];
else
    folder=source;
    requested_file=[];
end

% Create list of complete filenames (incl. path) in the folder
dirFolder = dir(folder);
filesInFolder = dirFolder(~[dirFolder.isdir]);
hidden = logical(ones(1,length(filesInFolder)));
for jj = 1:length(filesInFolder)
    if strcmp(filesInFolder(jj).name(1),'.')
        hidden(jj) = 0;
    end
end
filesInFolder = filesInFolder(hidden);%delete hidden files

if length(filesInFolder)==0
    % No files survived. This will not end well.
    error([ 'io_loadspec_dicom did not find any files in ' folder ])
end

requested_file_ix=find(strcmp(requested_file,{filesInFolder.name})); % NB, may be empty

if isfolder(folder)
    filesInFolder = fullfile(folder, {filesInFolder.name});
else
    % this can happen if the function is called with a wildcard, eg /path/to/abc*dcm
    filesInFolder = fullfile(dirname(folder), {filesInFolder.name});
end

% ARC170724 : If a single file is specified (rather than a directory), check
% whether this file should be considered in isolation (pre-averaged data) or
% if it should be combined with other items in the folder. We do this by
% inspecting the Series Number header.

if ~isempty(requested_file_ix)
    % check header of the first, last and requested items; if Series Number
    % does not match, don't try to combine
    check_header_ix = unique([1, requested_file_ix, length(filesInFolder) ]);
    encountered_series_numbers = [];
    encountered_problems = [];

    for ix = check_header_ix
        try
            dh = dicominfo(filesInFolder{ix});
            encountered_series_numbers(end+1) = dh.SeriesNumber;
        catch
            encountered_problems(end+1) = ix;
        end
    end

    encountered_series_numbers=unique(encountered_series_numbers);

    if any(encountered_problems == requested_file_ix)

        warning([ 'io_loadspec_dicom could not determine SeriesNumber of the requested file: ' source ]);
        % in this case, fall back on default behaviour

    elseif length(encountered_series_numbers)>1 || length(encountered_problems)>0

        % if we found more than one series, OR we found some unreadable files
        % (probably non-DICOM), then just retain that one specific file the
        % caller originally asked for

        if length(encountered_problems)>0
            warning([ 'io_loadspec_dicom encountered some non-DICOM data in ' source ]);
        end

        filesInFolder={source};
    end

    % in any other scenario, full back on the default behaviour
end

% Get the header of the first file to make some decisions.
DicomHeader = read_dcm_header(filesInFolder{1});
seqtype = DicomHeader.seqtype;

if isfield(DicomHeader, 'seqorig')
    seqorig = DicomHeader.seqorig;
else
% deal with missing seqorig field for dicom data load on Siemens Minnesota
% sequences (Auerbach and Deelchand versions, VB17A & VE11C)
    if contains(DicomHeader.sequenceFileName, 'slaser_dkd') %Deelchand/Oz
        seqorig = 'CMRR';
    elseif contains(DicomHeader.sequenceFileName, 'eja_svs') %Auerbach/Marjanska
        seqorig = 'CMRR';
    end
end


% Extract voxel dimensions
% If a parameter is set to zero (e.g. if no voxel rotation is
% performed), the respective field is left empty in the TWIX file. This
% case needs to be intercepted. Setting to the minimum possible value.
VoI_Params = {'VoI_InPlaneRot','VoI_RoFOV','VoI_PeFOV','VoIThickness','NormCor','NormSag','NormTra', ...
              'PosCor','PosSag','PosTra'};
for pp = 1:length(VoI_Params)
    if isempty(DicomHeader.(VoI_Params{pp}))
        DicomHeader.(VoI_Params{pp}) = realmin('double');
    end
end
geometry.size.VoI_RoFOV     = DicomHeader.VoI_RoFOV; % Voxel size in readout direction [mm]
geometry.size.VoI_PeFOV     = DicomHeader.VoI_PeFOV; % Voxel size in phase encoding direction [mm]
geometry.size.VoIThickness  = DicomHeader.VoIThickness; % Voxel size in slice selection direction [mm]
geometry.pos.PosCor         = DicomHeader.PosCor; % Coronal coordinate of voxel [mm]
geometry.pos.PosSag         = DicomHeader.PosSag; % Sagittal coordinate of voxel [mm]
geometry.pos.PosTra         = DicomHeader.PosTra; % Transversal coordinate of voxel [mm]
geometry.pos.TablePosTra    = DicomHeader.TablePosTra; % Transversal table position [mm]
geometry.rot.VoI_InPlaneRot = DicomHeader.VoI_InPlaneRot; % Voxel rotation in plane
geometry.rot.NormCor        = DicomHeader.NormCor; % Coronal component of normal vector of voxel
geometry.rot.NormSag        = DicomHeader.NormSag; % Sagittal component of normal vector of voxel
geometry.rot.NormTra        = DicomHeader.NormTra; % Transversal component of normal vector of voxel

% Preallocate array in which the FIDs are to be extracted.
if  DicomHeader.removeOS ~= 3
    fids = zeros(DicomHeader.vectorSize,length(filesInFolder));
else
    fids = zeros(DicomHeader.vectorSize*2,length(filesInFolder));
end
% Collect all FIDs and sort them into fids array
for kk = 1:length(filesInFolder)

    % First, attempt to retrieve the FID from the DICOM header:
    infoDicom = dicominfo(filesInFolder{kk});
    if isfield(infoDicom, 'SpectroscopyData')
        realFID = infoDicom.SpectroscopyData(1:2:end);
        imagFID = infoDicom.SpectroscopyData(2:2:end);
        fids(:,kk) = conj(realFID + 1j*imagFID);
    else
        % Try different route if that doesn't work:
        % Open DICOM
        fd = dicom_open(filesInFolder{kk});
        % read the signal in as a complex FID
        fids(:,kk) = dicom_get_spectrum_siemens(fd);
        fclose(fd);
    end

end


% Now make some decisions depending on sequence type, and the observed
% number of stored FIDs compared to the nominal number of acquisitions:
if contains(seqtype,'MEGA')
    % If the number of stored FIDs does not match twice the number of averages
    % stored in the DICOM header, the MEGA-PRESS data are averaged and
    % will (except for CMRR sequences) have the DIFF spectrum saved as well.
    % Delete this DIFF spectrum.
    nominalAvgs = 2*DicomHeader.nAverages;
    if size(fids,2) ~= nominalAvgs
        out.flags.averaged = 1;
        if ~strcmp(seqorig, 'CMRR')
            fids(:,3) = [];
        end
    else
        out.flags.averaged = 0;
    end

    % Currently, the DICOM recon of the universal sequence is flawed.
    % Kick out empty lines here and see if data can be reconstructed.
    if strcmp(seqorig, 'Universal') && out.flags.averaged == 0
        fids(:,size(fids,2)/2+1:end) = [];
    end

    % Rearrange into subspecs
    fids = reshape(fids,[size(fids,1) size(fids,2)/2 2]);

elseif strcmp(seqtype,'PRESS') || strcmp(seqtype,'STEAM') || strcmp(seqtype,'sLASER')
    % If the number of stored FIDs does not match the number of averages
    % stored in the DICOM header, the data are averaged.
    nominalAvgs = DicomHeader.nAverages;
    if (size(fids,2) ~= nominalAvgs) && (size(fids,2) == 1)
        out.flags.averaged = 1;
    else
        out.flags.averaged = 0;
    end

elseif strcmp(seqtype,'HERMES')

    out.flags.averaged = 0;
    % Currently, the DICOM recon of the universal sequence is flawed.
    % Kick out empty lines here and see if data can be reconstructed.
    if strcmp(seqorig, 'Universal') && out.flags.averaged == 0
        fids(:,size(fids,2)/2+1:end) = [];
    end
    fids_A = fids(:,1:4:end);
    fids_B = fids(:,2:4:end);
    fids_C = fids(:,3:4:end);
    fids_D = fids(:,4:4:end);
    fids = cat(3,fids_A,fids_B,fids_C,fids_D);
    dims.subSpecs=3;
end

% Assign correct dimensions
sz=size(fids);
Ncoils=1; % DICOM data are coil-combined on the scanner
Naverages = sz(2);

fids = (conj(fids));
sz=size(fids);
if ndims(fids)==4  %Default config when 4 dims are acquired
    dims.t=1;
    dims.coils=2;
    dims.averages=3;
    dims.subSpecs=4;
elseif ndims(fids)<4  %To many permutations...ask user for dims.
    if Naverages == 1 && Ncoils == 1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=3;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=0;
        end
    elseif Naverages>1 && Ncoils==1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=3;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=0;
        end
    elseif Naverages==1 && Ncoils>1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=3;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=0;
        end
    elseif Naverages>1 && Ncoils>1
        if ndims(fids)>3
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=4;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=0;
        end
    end
%     dims.t=1;
%     dims.coils=input('Enter the coils Dimension (0 for none):  ');
%     dims.averages=input('Enter the averages Dimension (0 for none):  ');
%     dims.subSpecs=input('Enter the subSpecs Dimension (0 for none);  ');
dims.extras = 0;
end

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired
%averages in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    if dims.averages~=0
        averages=sz(dims.averages)*sz(dims.subSpecs);
        rawAverages=averages;
    else
        averages=sz(dims.subSpecs);
        rawAverages=nominalAvgs;
    end
else
    if dims.averages~=0
        averages=sz(dims.averages);
        rawAverages=averages;
    else
        averages=1;
        rawAverages=nominalAvgs;
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

specs=fftshift(fft(fids,[],dims.t),dims.t);

%Now get relevant scan parameters:*****************************
Bo = DicomHeader.B0;
if  DicomHeader.removeOS ~= 3
    dwelltime = DicomHeader.dwellTime * 1e-9 * 2; % DICOM contain data with removed oversampling
else
    dwelltime = DicomHeader.dwellTime * 1e-9; % DICOM with oversampled data
end
%Calculate Dwell Time
spectralwidth=1/dwelltime;

%Get TxFrq
txfrq=DicomHeader.tx_freq;


%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=f/(Bo*42.577);

% Philips data assumes the center frequency to be 4.68 ppm:
centerFreq = 4.68;
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
out.seq=seqtype;
out.te=DicomHeader.TE;
out.tr=DicomHeader.TR;
out.pointsToLeftshift=0;
out.centerFreq = centerFreq;
out.geometry = geometry;
if isfield(DicomHeader, 'SoftwareVersions')
    out.software = DicomHeader.SoftwareVersions;
else
    out.software = 'Siemens';
end
%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.addedrcvrs=1;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end

% Add info for niiwrite
if isfield(DicomHeader,'PatientPosition')
    out.PatientPosition = DicomHeader.PatientPosition;
else
    out.PatientPosition = '';
end
out.Manufacturer = 'Siemens';
[~,filename,ext] = fileparts(filesInFolder{1});
out.OriginalFile = [filename ext];

% Sequence flags
out.flags.isUnEdited = 0;
out.flags.isMEGA = 0;
out.flags.isHERMES = 0;
out.flags.isHERCULES = 0;
out.flags.isPRIAM = 0;
out.flags.isMRSI = 0;
if strcmp(seqtype,'PRESS') || strcmp(seqtype,'STEAM') || strcmp(seqtype,'SLASER')
    out.flags.isUnEdited = 1;
end
if contains(seqtype,'MEGA')
    out.flags.isMEGA = 1;
end
if strcmp(seqtype,'HERMES')
    out.flags.isHERMES = 1;
end
if strcmp(seqtype,'HERCULES')
    out.flags.isHERCULES = 1;
end


%DONE
