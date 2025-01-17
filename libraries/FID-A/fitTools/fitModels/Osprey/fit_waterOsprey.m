% fit_waterOsprey.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [fitParams] = fit_waterOsprey(dataToFit, resBasisSet, fitOpts);
% 
% DESCRIPTION:
% This is the function describing the Osprey water fitting model, which is
% performed with the resampled water basis set resBasisSet on the data contained
% in the FID-A water data structure dataToFit, using the fit options in the
% structure fitOpts.
% 
% OUTPUTS:
% fitParams   = Structure containing all necessary fit parameters
%               x = concatenated vector with zero- and first-order phase
%                   shift, Gaussian linebroadening, Lorentzian
%                   linebroadening, and frequency shift
%               ampl = amplitude of the water basis function
%
% INPUTS:
% dataToFit   = FID-A data structure
% basisSet    = FID-A basis set container
% fitOpts     = Structure containing fit options

function [fitParams] = fit_waterOsprey(dataToFit, resBasisSet, fitOpts)


%%% 1. EXTRACT OPTIONS AND PREPARE FIT %%%
% Extract ppm fit range
fitRangePPMWater = fitOpts.rangeWater;


%%% 2. SET STARTING VALUES %%%
% Set the starting values for the non-linear parameters to be passed on
% to the Levenberg-Marquardt NLLS solving algorithm. Note that the
% amplitude parameter does not need to be initialized, as the
% non-negative linear Lawson-Hanson solver does not require starting
% values to be provided.
ph0         = 0; % zero-order phase correction [rad]
ph1         = 0; % first-order phase correction [rad]
gaussLB     = op_getLW(dataToFit,4.5,4.9); % Gaussian damping [Hz^2]
lorentzLB   = 0.01; % Lorentzian damping [Hz] for water
freqShift   = 0; % Frequency shift [Hz] for water
% Concatenate together into one large x0 vector.
x0          = [ph0; ph1; gaussLB; lorentzLB; freqShift];


%%% 3. CALL NON-LINEAR SOLVER %%%
% Run the non-linear solver to optimize the non-linear parameters. At each
% iteration of the non-linear solver, the linear parameters are calculated
% with a non-negative linear solver.
% Levenberg-Marquardt implementation allowing parameter constraints
% (c) Alexander Drentler
% https://www.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox
opts.Display            = 'notify';
opts.Jacobian           = 'limit';
%   lets get rid of all stopping criteria
opts.FooTol             = NaN; 
opts.RelFooTol          = NaN;
opts.FooTol             = NaN;
opts.RelTolX            = NaN;
opts.TolX               = 1e-10;
%   note that we need to square the acceptance tolerance
abstol                  = 1e-8;
opts.AccTol             = abstol^2;
opts.MaxDamping         = 1e20;
opts.Broyden_updates    = 'off';

if ~isfield(fitOpts, 'isMRSI')
    % Set constraints
    lb      = zeros(size(x0));
    ub      = zeros(size(x0));
    lb(1)   = -2*pi; 
    ub(1)   = 2*pi; % Zero order phase shift [rad]
    lb(2)   = -pi/4; 
    ub(2)   = pi/4; % First order phase shift [rad]
    lb(3)   = 0; 
    ub(3)   = 5000; % Gaussian linebroadening [Hz^2]
    lb(4)   = 0; 
    ub(4)   = 50; % Lorentzian linebroadening [Hz]
    lb(5)   = -15; % Frequency shift [Hz]
    ub(5)   = 15; % Frequency shift [Hz]
else
    % Set constraints
    lb      = zeros(size(x0));
    ub      = zeros(size(x0));
    lb(1)   = -2*pi; 
    ub(1)   = 2*pi; % Zero order phase shift [rad]
    lb(2)   = -pi; 
    ub(2)   = pi; % First order phase shift [rad]
    lb(3)   = 0; 
    ub(3)   = 5000; % Gaussian linebroadening [Hz^2]
    lb(4)   = 0; 
    ub(4)   = 50; % Lorentzian linebroadening [Hz]
    lb(5)   = -15; % Frequency shift [Hz]
    ub(5)   = 15; % Frequency shift [Hz]   
end

% Run non-linear solver
[x] = LM(@(x) fit_waterOsprey_Model(x, dataToFit, resBasisSet, fitRangePPMWater),x0,lb,ub,opts);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameter (i.e. the water amplitude).
% Apply the parameters from the final iteration of the non-linear solver to
% the basis set to create the phased&broadened water basis function.
[ampl, x] = fit_waterOsprey_FinalLinear(dataToFit, resBasisSet, x, fitRangePPMWater);


% Store all fit parameters in a structure
fitParams.ampl      = ampl;
fitParams.ph0       = x(1);
fitParams.ph1       = x(2);
fitParams.gaussLB   = x(3);
fitParams.lorentzLB = x(4);
fitParams.freqShift = x(5);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = fit_waterOsprey_Model(x, dataToFit, basisSet, fitRangePPMWater)
%% F = fit_nlinwaterOsprey(x, dataToFit, basisSet, fitRangePPMWater)
%   This function initializes and runs fitting of an unsuppressed water
%   signal with a simulated water basis function. Parameters that are
%   applied to the basis function are zero-order phase shift, frequency
%   shift, two parameters for Gaussian/Lorentzian lineshape contribution,
%   and amplitude.
%
%   The minimization problem is solved with a Levenberg-Marquardt
%   non-linear least-squares (NLLS) optimization for the nonlinear 
%   parameters (all except amplitude), with a fast non-negative linear 
%   least-squares (FNNLS) algorithm used at every iteration of the NLLS.
%
%   USAGE:
%       F           = fit_waterOsprey_Model(x, dataToFit, basisSet, fitRangePPMWater)
%
%   INPUTS:
%       data        = Vector containing the concatenated real and imaginary
%                     values of an acquired water signal.
%       basis       = Vector containing the concatenated real and imaginary
%                     values of a simulated water basis function.
%       x0          = Vector providing initial values to the non-linear
%                     least-squares solver.
%                     x0(1) = zero-order phase correction
%                     x0(2) = frequency shift
%                     x0(3) = parameter for Gaussian/Lorentzian lineshape
%                     x0(4) = parameter for Gaussian/Lorentzian lineshape
%
%   OUTPUTS:
%       F           = Objective function
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-09)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code uses the 'fnnls' implementation of a fast non-negative
%       least-squares optimization algorithm by Rasmus Bro.
%       https://www.mathworks.com/matlabcentral/fileexchange/3388-nnls-and-constrained-regression
%       Bro et al., Journal of Chemometrics 11:393-401 (1997)
%
%   HISTORY:
%       2019-04-09: First version of the code.

% 1. Cut out the frequency range of the spectrum to be fit
dataToFit = op_freqrange(dataToFit, fitRangePPMWater(1), fitRangePPMWater(end));
% Turn complex-valued problem into real-valued problem
data = [real(dataToFit.specs); imag(dataToFit.specs)];

% 2. Extract parameters
nBasisFcts  = 1;
ph0         = x(1); % zero-order phase correction [rad]
ph1         = x(2); % first-order phase correction [rad]
gaussLB     = x(3); % Gaussian damping [Hz^2]
lorentzLB   = x(4); % Lorentzian damping [Hz] for each basis function
freqShift   = x(5); % Frequency shift [Hz] for each basis function

% 3. Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(1i*freqShift.*t)' .* exp(-lorentzLB.*t)' .* exp(-gaussLB.*t.*t)';    
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% 4. Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph0) .* exp(1i*ph1*2*pi.*f)';    
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% Cut out the frequency range of the basis set
basisSet = op_freqrange(basisSet,fitRangePPMWater(1),fitRangePPMWater(end));


% 6. Run the first iteration of the non-negative linear solver
A = [real(basisSet.specs); imag(basisSet.specs)];
b = data;
[ampl]  = fnnls(A'*A,A'*b);

% 8. The functional to be used by the non-linear least squares solver is
% the fit residual:
F = b - (A*ampl);

end 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FINAL LINEAR ITERATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ampl, x] = fit_waterOsprey_FinalLinear(dataToFit, basisSet, x, fitRangePPMWater)

% 1. Extract parameters
ph0         = x(1); % zero-order phase correction [rad]
ph1         = x(2); % first-order phase correction [rad]
gaussLB     = x(3); % Gaussian damping [Hz^2]
lorentzLB   = x(4); % Lorentzian damping [Hz]
freqShift   = x(5); % Frequency shift [Hz]

% 2. Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
basisSet.fids = basisSet.fids .* exp(1i*freqShift.*t)' .* exp(-lorentzLB.*t)' .* exp(-gaussLB.*t.*t)';
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% 3. Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
basisSet.specs = basisSet.specs .* exp(1i*ph0) .* exp(1i*ph1*2*pi.*f)';
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% Cut out the frequency range of the basis set
appliedBasisSet = op_freqrange(basisSet, fitRangePPMWater(1), fitRangePPMWater(end));

% Run the final iteration of the non-negative linear solver
% Cut out the frequency range of the spectrum to be fit
dataToFit = op_freqrange(dataToFit, fitRangePPMWater(1), fitRangePPMWater(end));
% Turn complex-valued problem into real-valued problem
data = [real(dataToFit.specs); imag(dataToFit.specs)];
A = [real(appliedBasisSet.specs); imag(appliedBasisSet.specs)];
b = data;
[ampl]  = fnnls(A'*A,A'*b);

end 

