function out = osp_plotLoad(MRSCont, kk, which, stag, ppmmin, ppmmax, xlab, ylab, figTitle)
%% out = osp_plotLoad(MRSCont, kk, which, stag, ppmmin, ppmmax, xlab, ylab, figTitle)
%   Creates a figure showing raw data stored in an Osprey data container,
%   ie in the raw fields. This function will display the *unprocessed*
%   data, i.e. all averages will be shown prior to spectral alignment,
%   averaging, water removal, and other processing steps carried out in
%   OspreyProcess.
%
%   USAGE:
%       out = osp_plotLoad(MRSCont, kk, which, stag, ppmmin, ppmmax, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset
%       which    = String for the spectrum to fit
%                   OPTIONS:    'mets'
%                               'ref'
%                               'w'
%       stag     = Numerical value representing the fraction of the maximum
%                   by which individual averages should be plotted
%                   vertically staggered (optional. Default - 0, i.e. no stagger)
%       xlab     = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab     = label for the y-axis (optional.  Default = '');
%       figTitle = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    error('Trying to plot raw data, but no data has been loaded yet. Run OspreyLoad first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin<9
    switch which
        case 'mets'
            figTitle = sprintf('Load metabolite data plot');
        case 'ref'
            figTitle = sprintf('Load water reference data plot');
        case 'w'
            figTitle = sprintf('Load water data plot');
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
    if nargin<8
        ylab='';
        if nargin<7
            xlab='Frequency (ppm)';
            if nargin<6
                switch which
                    case 'mets'
                        ppmmax = 4.5;
                    case {'ref', 'w'}
                        ppmmax = 2*4.68;
                    otherwise
                        error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                end
                if nargin<5
                    switch which
                        case 'mets'
                            ppmmin = 0.2;
                        case {'ref', 'w'}
                            ppmmin = 0;
                        otherwise
                            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                    end
                    if nargin<4
                        stag = 0;
                        if nargin < 3
                            which = 'mets';
                            if nargin < 2
                                kk = 1;
                                if nargin<1
                                    error('ERROR: no input Osprey container specified.  Aborting!!');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw spectra in the plot range
switch which
    case 'mets'
        dataToPlot  = op_freqrange(MRSCont.raw{kk}, ppmmin, ppmmax);
    case 'ref'
        dataToPlot  = op_freqrange(MRSCont.raw_ref{kk}, ppmmin, ppmmax);
    case 'w'
        dataToPlot  = op_freqrange(MRSCont.raw_w{kk}, ppmmin, ppmmax);
    otherwise
        error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
end


%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
out = figure;
nAvgs = dataToPlot.averages;

% Add the data and plot
hold on;    
% Staggered plots will be in all black and separated by the mean of the
% maximum across all averages divided by the number of averages
stag = stag*max(abs(mean(max(real(dataToPlot.specs)))), abs(mean(min(real(dataToPlot.specs)))));
% Loop over all averages
for rr = 1:nAvgs
    plot(dataToPlot.ppm, dataToPlot.specs(:,rr) + rr*stag, 'k', 'LineWidth', 0.5);
end
hold off;


%%% 4. DESIGN FINETUNING %%%
set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
    set(gca, 'YColor', 'w');
else
    set(gca, 'YColor', 'k');
end
% Black axes, white background
set(gca, 'XColor', 'k');
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');
box off;
title(figTitle);
xlabel(xlab, 'FontSize', 16);
ylabel(ylab, 'FontSize', 16);


%%% 5. ADD OSPREY LOGO %%%
[I, map] = imread('osprey.gif','gif');
axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
imshow(I, map);
axis off;


end

   