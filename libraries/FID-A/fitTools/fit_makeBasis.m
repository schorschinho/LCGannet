% fit_makeBasis.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [BASIS] = fit_makeBasis(folder, sequence, addMMFlag)
% 
% DESCRIPTION:
% Generates a basis set in FID-A structure. The code will search all *.mat
% files in the input folder for FID-A structures with simulated spectra. It
% also performs sanity checks on the simulation parameters, and returns
% warnings if parameters are not identical for all parameters.
% 
% INPUTS:
% folder    = folder containing *.mat files representing FID-A structures
% sequence  = sequence type
%             OPTIONS:  'unedited' (default)
%                       'MEGA'
%                       'HERMES'
%                       'HERCULES'
% addMMFlag = Flag to decide whether MM and lipid basis functions should be
%               added to the basis set.
%             OPTIONS:  1 = Add MM+lip (Default)
%                       0 = Don't add MM+lip
%
% OUTPUTS:
% BASIS     = Simulated basis set in FID-A structure format. 

function [BASIS] = fit_makeBasis(folder, sequence, addMMFlag)

% Parse input arguments
if nargin < 3
    addMMFlag = 1;
    if nargin < 2
        sequence = 'unedited';
    end
end

% Collect *.mat filenames from input folder
mat_files       = dir([folder filesep '*.mat']);
mat_filenames   = strcat(folder, filesep, {mat_files.name});
idx = contains(mat_filenames, 'Ref');
mat_filenames(idx) = [];
nMets           = length(mat_filenames);

% Loop over all *.mat filenames, load their data, store in a buffer
for kk = 1:nMets
    temp = load(mat_filenames{kk});
    % Multiplexed experiments (e.g. MEGA/HERMES) have more than one sub-basis 
    % simulated per metabolite. Find out how many:
    basisFct = fieldnames(temp);
    % Load the signals, DC-correct and store them in separate dimensions
    for ll = 1:length(basisFct)
        temp.(basisFct{ll})            = op_dccorr(temp.(basisFct{ll}),'p');
        buffer.fids(:,kk,ll)           = temp.(basisFct{ll}).fids;
        buffer.specs(:,kk,ll)          = temp.(basisFct{ll}).specs;
    end
    % The following field should always be the same for all sub-bases
    % (unless something is seriously flawed with the simulation code)
    buffer.t(:,kk)              = temp.(basisFct{1}).t;
    buffer.ppm(:,kk)            = temp.(basisFct{1}).ppm;
    buffer.spectralwidth(kk)    = temp.(basisFct{1}).spectralwidth;
    buffer.dwelltime(kk)        = temp.(basisFct{1}).dwelltime;
    buffer.n(kk)                = temp.(basisFct{1}).sz(1);
    buffer.linewidth(kk)        = temp.(basisFct{1}).linewidth;
    buffer.Bo(kk)               = temp.(basisFct{1}).Bo;
    buffer.seq{kk}              = temp.(basisFct{1}).seq;
    if isfield(temp.(basisFct{1}),'name')
        buffer.name{kk}             = temp.(basisFct{1}).name;
    else
        C = strsplit(mat_files(kk).name,'_');
        C = C{end};
        buffer.name{kk} = strrep(C,'.mat','');
    end
    buffer.te(kk)               = temp.(basisFct{1}).te;
    buffer.dims                 = temp.(basisFct{1}).dims;
    buffer.flags                = temp.(basisFct{1}).flags;
    if isfield(temp.(basisFct{1}), 'centerFreq')
        buffer.centerFreq           = temp.(basisFct{1}).centerFreq;
    else
        buffer.centerFreq = 3;
    end
end

% If the ppm axis runs from negative to positive, this has run with the old
% FID-A simulations - correct here:
% if issorted(buffer.ppm, 'descend')
%         sw = buffer.spectralwidth;
%         points = length(buffer.t);
%         freq = [(-sw/2)+(sw/(2*points)):sw/(points):(sw/2)-(sw/(2*points))];
%         ppm = freq/(buffer.Bo(1)*42.577);
%         ppm = ppm - 4.68;
%         % invert ppm axis
%         buffer.ppm = ppm;
%         % re-scale and fft FIDs
%        buffer.ppm = buffer.ppm + 0.04;
        buffer.specs = fftshift(fft(buffer.fids,[],1),1);
%end

% Test whether parameters are the same across all basis functions; flag
% warning if they are not; write into basis set struct if they are.
seq_params = {'spectralwidth','dwelltime','n','linewidth','Bo','seq','te', 'centerFreq'};
for pp = 1:length(seq_params)
    unique_params = unique(buffer.(seq_params{pp}));
    if length(unique_params) > 1
        error('WARNING! One or more sequence parameters are not the same across all input basis functions.');
    else
        BASIS.(seq_params{pp}) = unique_params;
    end
end

% Test whether ppm and t aves are the same across all basis functions; flag
% error if they are not; write into basis set struct if they are.
seq_params = {'ppm','t'};
for pp = 1:length(seq_params)
    unique_params = unique(buffer.(seq_params{pp}),'stable');
    if length(unique_params) ~= BASIS.n
        error('WARNING! One or more sequence parameters are not the same across all input basis functions.');
    else
        BASIS.(seq_params{pp}) = unique_params';
    end
end

% If chosen, add MM
if addMMFlag
    n = BASIS.n;
    sw = BASIS.spectralwidth;
    Bo = BASIS.Bo;
    lw = BASIS.linewidth;
    centerFreq = BASIS.centerFreq;
    % The amplitude and FWHM values are determined as for the LCModel and
    % TARQUIN algorithms (see Wilson et al., MRM 2011).
    hzppm = Bo*42.577;
    % op_gaussianPeak with amp = 1 produces a signal with the amplitude
    % equivalent to two fully simulated protons
    MM09            = op_gaussianPeak(n,sw,Bo,centerFreq,0.14*hzppm+lw,0.91,3/2);
    MMBase.MM09     = op_dccorr(MM09,'p');
    MM12            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm+lw,1.21,2/2);
    MMBase.MM12     = op_dccorr(MM12,'p');
    MM14            = op_gaussianPeak(n,sw,Bo,centerFreq,0.17*hzppm+lw,1.43,2/2);
    MMBase.MM14     = op_dccorr(MM14,'p');
    MM17            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm+lw,1.67,2/2);
    MMBase.MM17     = op_dccorr(MM17,'p');
    MM20a           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm+lw,2.08,1.33/2);
    MM20b           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm+lw,2.25,0.33/2);
    MM20c           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm+lw,1.95,0.33/2);
    MM20d           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm+lw,3.0,0.4/2);
    MM20            = op_addScans(MM20a,MM20b); MM20 = op_addScans(MM20,MM20c); MM20 = op_addScans(MM20,MM20d);
    MMBase.MM20     = op_dccorr(MM20,'p');
    Lip09           = op_gaussianPeak(n,sw,Bo,centerFreq,0.14*hzppm+lw,0.89,3/2);
    MMBase.Lip09    = op_dccorr(Lip09,'p');
    Lip13a          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm+lw,1.28,2/2);
    Lip13b          = op_gaussianPeak(n,sw,Bo,centerFreq,0.89*hzppm+lw,1.28,2/2);
    Lip13           = op_addScans(Lip13a,Lip13b);
    MMBase.Lip13    = op_dccorr(Lip13,'p');
    Lip20a          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm+lw,2.04,1.33/2);
    Lip20b          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm+lw,2.25,0.33/2);
    Lip20c          = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm+lw,2.8,0.33/2);
    Lip20           = op_addScans(Lip20a,Lip20b); Lip20 = op_addScans(Lip20,Lip20c);
    MMBase.Lip20    = op_dccorr(Lip20,'p');
    MMLips = {'MM09','MM12','MM14','MM17','MM20','Lip09','Lip13','Lip20'};
    
    % Now copy over the names, fids, and specs into the basis set structure
    for rr = 1:length(MMLips)
        buffer.name{nMets+rr}       = MMLips{rr};
        for qq = 1:size(length(basisFct))
            buffer.fids(:,nMets+rr,qq)   = MMBase.(MMLips{rr}).fids;
            buffer.specs(:,nMets+rr,qq)  = MMBase.(MMLips{rr}).specs;
        end
    end
    
    BASIS.flags.addedMM     = 1;
    BASIS.nMM               = length(MMLips);
    save_str = '_MM';
else
    BASIS.flags.addedMM     = 0;
    BASIS.nMM               = 0;
    save_str = '_noMM';
end

% If spectral editing has been performed, do the SUM and DIFF spectra here
% TO DO: CREATE AUTOMATIC RECOGNITION OF ON/OFF
if strcmp(sequence, 'MEGA')
    buffer.fids(:,:,3)      = buffer.fids(:,:,2) - buffer.fids(:,:,1); % DIFF
    buffer.specs(:,:,3)     = buffer.specs(:,:,2) - buffer.specs(:,:,1);
    buffer.fids(:,:,4)      = buffer.fids(:,:,2) + buffer.fids(:,:,1); % SUM
    buffer.specs(:,:,4)     = buffer.specs(:,:,2) + buffer.specs(:,:,1);
elseif strcmp(sequence, 'HERMES') || strcmp(sequence, 'HERCULES')
    buffer.fids(:,:,5)      = buffer.fids(:,:,2) + buffer.fids(:,:,3) - buffer.fids(:,:,1) - buffer.fids(:,:,4); % DIFF1 (GABA)
    buffer.specs(:,:,5)     = buffer.specs(:,:,2) + buffer.specs(:,:,3) - buffer.specs(:,:,1) - buffer.specs(:,:,4);
    buffer.fids(:,:,6)      = buffer.fids(:,:,1) + buffer.fids(:,:,3) - buffer.fids(:,:,2) - buffer.fids(:,:,4); % DIFF2 (GSH)
    buffer.specs(:,:,6)     = buffer.specs(:,:,1) + buffer.specs(:,:,3) - buffer.specs(:,:,2) - buffer.specs(:,:,4);
    buffer.fids(:,:,7)      = buffer.fids(:,:,1) + buffer.fids(:,:,3) + buffer.fids(:,:,2) + buffer.fids(:,:,4); % SUM
    buffer.specs(:,:,7)     = buffer.specs(:,:,1) + buffer.specs(:,:,3) + buffer.specs(:,:,2) + buffer.specs(:,:,4);
end

% Copy over the FID, specs, dims, and the metabolite names
BASIS.fids              = buffer.fids;
BASIS.specs             = buffer.specs;
BASIS.name              = buffer.name;
BASIS.dims              = buffer.dims;
BASIS.flags             = buffer.flags;
BASIS.nMets             = nMets;
BASIS.sz                = size(BASIS.fids);

% Save as *.mat file
save(['BASIS' save_str '.mat'], 'BASIS');

end

   