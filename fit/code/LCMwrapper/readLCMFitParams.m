function fitParams = readLCMFitParams(MRSCont, which, kk)

% Fall back to defaults if not provided
if nargin < 3
    which = 'A';
    if nargin < 2
        kk = 1;
        if nargin < 1
            error('ERROR: no input Osprey container specified.  Aborting!!');
        end
    end
end

% Determine output file strings
switch which
    case {'A','B','C','D'}
        lcmOutputFile = ['outputfile' which];
    case {'diff1'}
        lcmOutputFile = 'outputfileDiff1';
    case {'diff2'}
        lcmOutputFile = 'outputfileDiff2';
    case {'sum'}
        lcmOutputFile = 'outputfileSum';
end

% Read the fit results from the .table files
if isfile([MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk} '.table'])
    tab                 = mrs_readLcmodelTABLE([MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk} '.table']);
else
    error('ERROR! Cannot find LCModel output table file %s. \nVery likely, LCModel has not been executed correctly.\n', [MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk} '.table'])
end
fitParams.name      = tab.name;
fitParams.CRLB      = tab.SDpct;
fitParams.relConc   = tab.relative_conc;
fitParams.ph0       = tab.ph0;
fitParams.ph1       = tab.ph1;
fitParams.refShift  = tab.refShift;
fitParams.refFWHM   = tab.fwhm;
fitParams.SNR       = tab.snr;

%Remove the - in -CrCH2 because it interferes with the downstream functions
idx = find(strcmp(fitParams.name,'-CrCH2'));
if ~isempty(idx)
    fitParams.name{idx} = 'CrCH2';
end
%Remove the + in combinations because it interferes with the downstream functions
idx = find(contains(fitParams.name,'+'));
if ~isempty(idx)
    for combs = 1 : length(idx)
        fitParams.name{idx(combs)} = strrep(fitParams.name{idx(combs)},'+','_');
    end
end

% Read the spectrum, fit, and baseline from the .coord files
[ spectra, spectra_metabolites, x_ppm, info ] = mrs_readLcmodelCOORD( [MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk} '.coord'] );
fitParams.ppm           = x_ppm;
fitParams.data          = spectra(:,1);
fitParams.completeFit   = spectra(:,2);
if size(spectra,2) == 3
    fitParams.baseline      = spectra(:,3);
else
    fitParams.baseline      =zeros(size(spectra,1),1);
end
fitParams.residual      = fitParams.data - fitParams.completeFit;

% Add Nan values for nicer plots
if ~isempty(MRSCont.opts.fit.GAP.(which))
    idx = find(x_ppm < MRSCont.opts.fit.GAP.(which)(1));
    idx_1 = idx(1);
    idx_2 = idx_1-1;
    fitParams.data(idx_1) = NaN;
    fitParams.data(idx_2) = NaN;
    fitParams.completeFit(idx_1) = NaN;
    fitParams.completeFit(idx_2) = NaN;
    fitParams.baseline(idx_1) = NaN;
    fitParams.baseline(idx_2) = NaN;
    fitParams.residual(idx_1) = NaN;
    fitParams.residual(idx_2) = NaN;
    spectra_metabolites(idx_1,:) = NaN;
    spectra_metabolites(idx_2,:) = NaN;
end

% The .coord files also contain the individual metabolite fits, BUT only if
% the estimate is not zero, and the individual metabolite fits include the
% baseline.
fitParams.indivMets     = zeros(info.n, length(fitParams.name));
for rr = 1:length(fitParams.name)
    % Check whether a particular metabolite has been fit
    idxMatch = find(strcmp(info.metabolites, fitParams.name{rr}));
    % If it has been fit, save the fit with the correct index after
    % subtracting the baseline
    if idxMatch
        fitParams.indivMets(:,rr) = spectra_metabolites(:,idxMatch) - fitParams.baseline;
    else
        if ~isempty(MRSCont.opts.fit.GAP.(which))
            fitParams.indivMets(idx_1,rr) = NaN;
            fitParams.indivMets(idx_2,rr) = NaN;
        end
    end
end
% Store amplitudes regardless whether the are fit or not
fitParams.ampl        = tab.concentration';

% Read the raw area of the unsuppressed water peak from the .print file
infoPrint = mrs_readLcmodelPRINT( [MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk} '.print'] );
if isfield(infoPrint,'h2oarea')
    fitParams.h2oarea = infoPrint.h2oarea;
end



end