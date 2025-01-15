function [NoiseSD]=fit_getNoiseSD(in,noiseppmmin,noiseppmmax);
% [NoiseSD]=osp_gLCM_getNoiseSD(in,noiseppmmin,noiseppmmax);
% Calcualte the standard deviation of the noise using the noise covariance
% matrix
%
% USAGE:
% [NoiseSD]=osp_gLCM_getNoiseSD(in,noiseppmmin,noiseppmmax);
%
%   INPUTS:
%       in             = specs vector
%       noiseppmmin    = min of frequency range in which to measure noise.
%                       (Optional.  Default = 8.5 ppm);
%       noiseppmmax    = max of frequency range in which to measure noise.
%                       (Optional.  Default = 12 ppm);
%
%   OUTPUTS:
%       Noise            = Estimated standard deviation of noise.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu

%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017) 
%%  Diverge to default options if required 
if nargin<3
    noiseppmmax=12;                                                         % Set ppm range max
    if nargin<2
        noiseppmmin=8.5;                                                    % Set ppm range min
    end
end
%% Fail save if the frequency axis is too small

if max(in.ppm) < noiseppmmax                                                % Noise window is outside of data ppm window
    noiseppmmax = max(in.ppm)-0.1;                                          % Find maximum possible ppm minus 0.1 ppm to avoid edge effects
    noiseppmmin = noiseppmmax-3.5;                                          % Set new minimum value to max minus 3.5 to be similar to default length
end
%% Calculate noise estimates
ppmwindow=in.ppm(in.ppm>noiseppmmin & in.ppm<noiseppmmax)';                 % Exctract frequency range
noisewindow=in.specs(in.ppm>noiseppmmin & in.ppm<noiseppmmax,:);            % Extract noise range

if ~in.dims.extras                                                          % 1D data
    P=polyfit(ppmwindow,noisewindow,2);                                     % Calculate second order polynom
    noise=noisewindow-polyval(P,ppmwindow);                                 % Detrend noise
else
    noise = zeros(length(noisewindow),in.sz(in.dims.extras));               % Setup noise vector 
    for kk = 1 : in.sz(in.dims.extras)                                      % Loop over indirect dimension   
        P=polyfit(ppmwindow,noisewindow(:,kk),2);                           % Calculate second order polynom
        noise(:,kk)=noisewindow(:,kk)-polyval(P,ppmwindow);                 % Detrend noise
    end                                                                     % End loop over indirect dimension
end
noisecovariance = cov(real(noise));                                         % Calculate covariance matrix
NoiseSD = sqrt(diag(noisecovariance))';                                      % Calculate standard deviation
end



