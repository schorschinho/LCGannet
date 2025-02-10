function DicomHeader = read_dcm_header(fid)
%% dcmHeader = read_dcm_header(fid)
%   Reads header information from a DICOM file.
%
%   Example:
%       dcmHeader = read_dcm_header('file0001.dcm')
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-04-24)
%       goeltzs1@jhmi.edu
%   
%   Credits:    
% 
%   Version history:
%   0.9:  First version (2018-04-24)
%   0.91: Several sequence-specific loading fixes (2018-05-13)
%   0.92: Added support for sLASER sequence (2018-07-18)
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HEADER INFO PARSING %%%
% Simply open the dicom file, the information should all be the same.
fid = fopen(fid);

% Start looking for a convenient parameter block. The line before will 
% start with ### ASCCONV BEGIN and end with ### ASCCONV END. 
% This is defined here.
head_start_text = '### ASCCONV BEGIN';
head_end_text   = '### ASCCONV END';
tline = fgets(fid); % get first line



% Keep looking until start of the parameter block is found.
while (isempty(strfind(tline, head_start_text)))
    tline = fgets(fid);
    if contains(tline,'HFS')
        dcmHeader.PatientPosition = 'HFS';
    end
    if contains(tline,'syngo')
        ind = strfind(tline,'syngo');
        ind2 = strfind(tline,'LO');
        border = find(ind2 > ind);
        try
            dcmHeader.SoftwareVersions = tline(ind : ind2(border(1))-2);
        end
    end
end

% Look for regular expression containing the 'equal' signs
while (isempty(strfind(tline, head_end_text)))
    [tokens,matches] = regexp(tline,'([\w\[\].]*)\s*=\s*([\w.-\"\\]*)','tokens','match');
    % When a matching string is found, parse the results into a struct
    if length(tokens) == 1
        fieldname = regexprep(tokens{1}{1}, '\[|\]|_',''); % delete invalid characters
        if isempty(strfind(tokens{1}{2},'"'))
            if strcmp(tokens{1}{2},'0x1')
                value = 0;
            elseif strcmp(tokens{1}{2},'1x1')
                value = 1;
            else
                value = str2double(tokens{1}{2});
            end
        else
            value = tokens{1}{2};
        end
        C = strsplit(fieldname,'.'); % check for nested variable names
        switch length(C)
            case 1
                % Intercept and skip read errors
                try
                    dcmHeader.(C{1}) = value;
                catch
                    tline = fgets(fid);
                    continue;
                end
            case 2
                dcmHeader.(C{1}).(C{2}) = value;
            case 3
                dcmHeader.(C{1}).(C{2}).(C{3}) = value;
            case 4
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}) = value;
            case 5
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}).(C{5}) = value;
            case 6
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}).(C{5}).(C{6}) = value;
            case 7
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}).(C{5}).(C{6}).(C{7}) = value;   
        end
    end
    tline = fgets(fid);
end
fclose(fid);


% Determine sequence name and type
DicomHeader.sequenceFileName     = dcmHeader.tSequenceFileName; % Full sequence name
% Determine the origin of the sequence
if contains(DicomHeader.sequenceFileName,'svs_edit')
    DicomHeader.seqtype = 'MEGA_PRESS';
    if contains(DicomHeader.sequenceFileName,'univ') %#ok<STRIFCND>
        DicomHeader.seqorig = 'Universal'; % Universal sequence
        if (dcmHeader.sWipMemBlock.alFree7 == 1)
            DicomHeader.seqtype = 'HERMES';
        end
    else
        DicomHeader.seqorig = 'WIP'; % Siemens WIP
    end
elseif contains(DicomHeader.sequenceFileName,'jn_')
    DicomHeader.seqtype = 'MEGA';
    DicomHeader.seqorig = 'JN'; % Jamie Near's sequence
elseif (contains(DicomHeader.sequenceFileName,'eja_svs_mpress') || contains(DicomHeader.sequenceFileName,'eja_svs_mslaser'))
   if contains(DicomHeader.sequenceFileName,'eja_svs_mpress')
        DicomHeader.seqtype = 'MEGA_PRESS';
   else
        DicomHeader.seqtype = 'MEGA_SLASER';
   end
    DicomHeader.seqorig = 'CMRR'; % Minnesota sequence
elseif contains(DicomHeader.sequenceFileName,'svs_se')
    DicomHeader.seqtype = 'PRESS'; % PRESS
    DicomHeader.seqorig = 'unknown'; % Unknwon
elseif contains(DicomHeader.sequenceFileName,'eja_svs')
    DicomHeader.seqtype = 'PRESS'; % PRESS
    DicomHeader.seqorig = 'CMRR'; % Minnesota sequence
elseif contains(DicomHeader.sequenceFileName,'svs_slaser')
    DicomHeader.seqtype = 'sLASER'; % sLASER
    DicomHeader.seqorig = 'unknown'; % Unknwon
elseif contains(DicomHeader.sequenceFileName,'svs_st')
    DicomHeader.seqtype = 'STEAM'; % STEAM
    if contains(DicomHeader.sequenceFileName,'vapor')
       DicomHeader.seqorig = 'CMRR'; % Probably CMRR
    else
       DicomHeader.seqorig = 'WIP'; % Probably WIP or product?
    end
else
    DicomHeader.seqorig = DicomHeader.sequenceFileName;
    error(['Unknown sequence: ' DicomHeader.seqorig '. Please consult the Gannet team for support.'])
end

% Read information
DicomHeader.TR                   = dcmHeader.alTR0 * 1e-3; % TR [ms]
DicomHeader.TE                   = dcmHeader.alTE0 * 1e-3; % TE [ms]
if isfield(dcmHeader, 'lAverages')
    DicomHeader.nAverages        = dcmHeader.lAverages;
else
    % Minnesota sequence (CMRR, Eddy Auerbach) may store numbers of averages in a
    % different field. GO 112017. Spelling may vary as well...
    if isfield(dcmHeader, 'sWipMemBlock')
        DicomHeader.nAverages        = dcmHeader.sWipMemBlock.alFree2;
    elseif isfield(dcmHeader, 'sWiPMemBlock')
        DicomHeader.nAverages        = dcmHeader.sWiPMemBlock.alFree2;
    end
end
try
    DicomHeader.removeOS             = dcmHeader.sSpecPara.ucRemoveOversampling; % Is the oversampling removed in the RDA files?
catch
    DicomHeader.removeOS             = 3; % We assume that oversampling is still applied
end
DicomHeader.vectorSize           = dcmHeader.sSpecPara.lVectorSize; % Data points specified on exam card
% GO180424: If a parameter is set to zero (e.g. if no voxel rotation is
% performed), the respective field does not show up in the dicom file. This
% case needs to be intercepted. Setting to the minimum possible value.
if ~isfield(dcmHeader.sSpecPara.sVoI, 'dInPlaneRot')
    dcmHeader.sSpecPara.sVoI.dInPlaneRot = realmin('double');
end
if ~isfield(dcmHeader.sSpecPara.sVoI, 'sPosition')
    dcmHeader.sSpecPara.sVoI.sPosition.dCor = realmin('double');
end
VoI_Params = {'dCor','dSag','dTra'};
for pp = 1:length(VoI_Params)
    if ~isfield(dcmHeader.sSpecPara.sVoI.sNormal, VoI_Params{pp})
        dcmHeader.sSpecPara.sVoI.sNormal.(VoI_Params{pp}) = realmin('double');
    end
    if ~isfield(dcmHeader.sSpecPara.sVoI.sPosition, VoI_Params{pp})
        dcmHeader.sSpecPara.sVoI.sPosition.(VoI_Params{pp}) = realmin('double');
    end
end

DicomHeader.VoI_InPlaneRot       = dcmHeader.sSpecPara.sVoI.dInPlaneRot; % Voxel rotation in plane
DicomHeader.VoI_RoFOV            = dcmHeader.sSpecPara.sVoI.dReadoutFOV; % Voxel size in readout direction [mm]
DicomHeader.VoI_PeFOV            = dcmHeader.sSpecPara.sVoI.dPhaseFOV; % Voxel size in phase encoding direction [mm]
DicomHeader.VoIThickness         = dcmHeader.sSpecPara.sVoI.dThickness; % Voxel size in slice selection direction [mm]
DicomHeader.NormCor              = dcmHeader.sSpecPara.sVoI.sNormal.dCor; % Coronal component of normal vector of voxel
DicomHeader.NormSag              = dcmHeader.sSpecPara.sVoI.sNormal.dSag; % Sagittal component of normal vector of voxel
DicomHeader.NormTra              = dcmHeader.sSpecPara.sVoI.sNormal.dTra; % Transversal component of normal vector of voxel
DicomHeader.PosCor               = dcmHeader.sSpecPara.sVoI.sPosition.dCor; % Coronal coordinate of voxel [mm]
DicomHeader.PosSag               = dcmHeader.sSpecPara.sVoI.sPosition.dSag; % Sagittal coordinate of voxel [mm]
DicomHeader.PosTra               = dcmHeader.sSpecPara.sVoI.sPosition.dTra; % Transversal coordinate of voxel [mm]
DicomHeader.TablePosTra          = dcmHeader.lScanRegionPosTra; % Transversal table position [mm]

% delta frequency (center of slice selection)
if isfield(dcmHeader.sSpecPara, 'dDeltaFrequency')
    DicomHeader.deltaFreq        = dcmHeader.sSpecPara.dDeltaFrequency;
else
    DicomHeader.deltaFreq        = 0;
end

DicomHeader.dwellTime            = dcmHeader.sRXSPEC.alDwellTime0; % dwell time [ns]
DicomHeader.tx_freq              = dcmHeader.sTXSPEC.asNucleusInfo0.lFrequency; % Transmitter frequency [Hz]

if isfield(dcmHeader.sProtConsistencyInfo, 'flNominalB0')
    DicomHeader.B0                   = dcmHeader.sProtConsistencyInfo.flNominalB0; % Nominal B0 [T]
else % XA30 missing this field, so determine B0 from Tx freq:
    Nucleus = dcmHeader.sTXSPEC.asNucleusInfo0.tNucleus;
    Nucleus = Nucleus(isstrprop(Nucleus,'alphanum'));

    switch strtrim(Nucleus) % Determine correct Gyromagnetic ratio for the Nucleus (3 s.f.):
        case '1H'
            gamma = 42.577;
        case '2H'
            gamma = 6.536;
        case '3HE'
            gamma = -32.434;
        case '7LI'
            gamma = 16.546;
        case '13C'
            gamma = 10.708;
        case '19F'
            gamma = 40.052;
        case '23NA'
            gamma = 11.262;
        case '31P'
            gamma = 17.235;
        case '129XE'
            gamma = -11.777;
    end
    DicomHeader.B0 = DicomHeader.tx_freq*1e-6 ./ gamma; % determine B0 from the transmit frequency [T]
end

if isfield(dcmHeader,'PatientPosition')
    DicomHeader.PatientPosition = dcmHeader.PatientPosition;
end
if isfield(dcmHeader,'SoftwareVersions')
    DicomHeader.SoftwareVersions = dcmHeader.SoftwareVersions;
end
% these may only be extractable from a few sequences and MEGA-PRESS
% versions:
% editing pulse parameters
if isfield(DicomHeader, 'seqorig')
    if strcmp(DicomHeader.seqorig,'CMRR')
        if contains(DicomHeader.seqtype, 'MEGA')
            if isfield(dcmHeader, 'sWipMemBlock')
                if isfield(dcmHeader.sWipMemBlock, 'adFree3')
                    DicomHeader.editRF.freq(1) = dcmHeader.sWipMemBlock.adFree3;
                end
                if isfield(dcmHeader.sWipMemBlock, 'adFree2')
                    DicomHeader.editRF.freq(2) = dcmHeader.sWipMemBlock.adFree2;
                end
                if isfield(dcmHeader.sWipMemBlock, 'adFree6')
                    DicomHeader.editRF.bw = dcmHeader.sWipMemBlock.adFree8;
                end
            elseif isfield(dcmHeader, 'sWiPMemBlock')
                if isfield(dcmHeader.sWiPMemBlock, 'adFree3')
                    DicomHeader.editRF.freq(1) = dcmHeader.sWiPMemBlock.adFree3;
                end
                if isfield(dcmHeader.sWiPMemBlock, 'adFree2')
                    DicomHeader.editRF.freq(2) = dcmHeader.sWiPMemBlock.adFree2;
                end
                if isfield(dcmHeader.sWiPMemBlock, 'adFree6')
                    DicomHeader.editRF.bw = dcmHeader.sWiPMemBlock.adFree8;
                end
            end
        else
            % For now, we probably don't need much from this block.
        end
    else
        if isfield(dcmHeader, 'sWipMemBlock')
            if isfield(dcmHeader.sWipMemBlock, 'adFree9')
                DicomHeader.editRF.centerFreq = dcmHeader.sWipMemBlock.adFree9;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree7')
                DicomHeader.editRF.freq(1) = dcmHeader.sWipMemBlock.adFree7;
                DicomHeader.editRF.freq(2) = DicomHeader.editRF.centerFreq + (DicomHeader.editRF.centerFreq - DicomHeader.editRF.freq(1));
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree8')
                DicomHeader.editRF.bw = dcmHeader.sWipMemBlock.adFree8;
            end
        elseif isfield(dcmHeader, 'sWiPMemBlock')
            if isfield(dcmHeader.sWiPMemBlock, 'adFree9')
                DicomHeader.editRF.centerFreq = dcmHeader.sWiPMemBlock.adFree9;
            end
            if isfield(dcmHeader.sWiPMemBlock, 'adFree7')
                DicomHeader.editRF.freq(1) = dcmHeader.sWiPMemBlock.adFree7;
                DicomHeader.editRF.freq(2) = DicomHeader.editRF.centerFreq + (DicomHeader.editRF.centerFreq - DicomHeader.editRF.freq(1));
            end
            if isfield(dcmHeader.sWiPMemBlock, 'adFree8')
                DicomHeader.editRF.bw = dcmHeader.sWiPMemBlock.adFree8;
            end
        end
    end
end

