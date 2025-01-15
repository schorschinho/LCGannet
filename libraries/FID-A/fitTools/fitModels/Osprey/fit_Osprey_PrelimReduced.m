function [fitParamsStep1] = fit_Osprey_PrelimReduced(dataToFit, basisSet, fitRangePPM)
%% [fitParamsStep1] = fit_Osprey_PrelimReduced(dataToFit, basisSet, fitRangePPM)
%   Performs the first step of the LCModel preliminary analysis
%   analogous to the LCModel algorithm. The algorithm is described in:
%       S.W. Provencher, "Estimation of metabolite concentrations from
%       localized in vivo NMR spectra", Magn Reson Med 30(6):672-679 (1993)
%
%   During the first step, the input spectrum is fit using a reduced basis
%   set (Cr, Glu, mI, GPC, NAA) and simplified model, using a common
%   frequency shift and common Gaussian and Lorentzian linebroadening
%   for all basis functions of the reduced basis set.
%
%   In addition, there is an unregularized baseline contribution. It seems
%   like this initial baseline has a fixed number of splines (18?), but I 
%   need to investigate this further. For now, set the minKnotSpacingPPM
%   parameter (= DKNTMN in LCModel) to 0.16.
%
%   Input:
%       dataToFit       = FID-A data structure
%       basisSet        = FID-A basis set container
%       fitRangePPM     = 2-element vector: fit range [ppm]
%                         (this is the range over which the difference 
%                         between spectrum and model is minimized)
%   Output:
%       fitParamsStep1  = Fit parameters:
%                         - amplitudes of basis functions
%                         - zero-order phase       [deg]
%                         - first-order phase      [deg/ppm]
%                         - Gaussian LB            [Hz]
%                         - Lorentzian LB          [Hz]
%                         - global frequency shift [Hz]
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2020-01-14)
%       goeltzs1@jhmi.edu
%
%   History:
%       2020-01-14: First version of the code.
%


%%% 0. SET EXPECTATION VALUES FOR LORENTZIAN LINEBROADENING AND FREQ SHIFTS
% See LCModel manual, Section 11.14 (p. 149)
exT2 = 2.0;     % [Hz] - expectation value for Lorentzian linebroadening
SDT2 = 0.4;     % [Hz] - standard deviation for Lorentzian linebroadening
SDSH = 0.004;   % [Hz] - standard deviation for frequency shifts
scalingT2 = sqrt(dataToFit.txfrq*1e-6 / 85.15); % scaling factor to account for T2 decrease with field strength


%%% 1. CREATE REDUCED BASIS SET %%%
% By default, the reduced basis set includes only NAA, Cr, PCh, mI, and
% Glu, analogous to LCModel. If those are not available creatine- and
% choline-containing basis functions are used. If mI and Glu are not
% included in the basis set theywill be omitted in the reduced fit.

if ~isempty(find(strcmp(basisSet.name, 'NAA')))
    metabList.NAA   = 1;
end
if ~isempty(find(strcmp(basisSet.name, 'Cr')))
    metabList.Cr    = 1;
else if ~isempty(find(strcmp(basisSet.name, 'PCr')))
    metabList.PCr    = 1;
    end
end

if ~isempty(find(strcmp(basisSet.name, 'PCh')))
    metabList.PCh    = 1;
else if ~isempty(find(strcmp(basisSet.name, 'GPC')))
    metabList.GPC    = 1;
    else if ~isempty(find(strcmp(basisSet.name, 'fCho')))
            metabList.fCho    = 1;
        end
    end
end

if ~isempty(find(strcmp(basisSet.name, 'mI')))
    metabList.mI   = 1;
end

if ~isempty(find(strcmp(basisSet.name, 'Glu')))
    metabList.Glu   = 1;
end


if ~exist('metabList','var')
    msg = 'Your metabolite list does not contain any of the standard metabolites (NAA, Cr, PCh, mI, or Glu) used for the preliminary fiting. You have to either include one of those or you have to update the prelinimary fitting function';
    fprintf(msg);
    error(msg);
end

nMets = numel(fieldnames(metabList));
fitMM = 0;
reducedBasisSet     = fit_selectMetabs(basisSet, metabList, fitMM);
% Create the spline basis functions for the given resolution, fit range,
% and knot spacing parameter.
[splineArray]       = fit_makeSplineBasis(dataToFit, fitRangePPM, 0.16);
nSplines            = size(splineArray,2);


%%% 2. SET AND GET STARTING VALUES %%%
% Set the starting values for all parameters.
ph0  = 0; % zero-order phase correction [deg]
ph1  = 0; % first-order phase correction [deg/ppm]
gaussLB     = 0.04 * dataToFit.txfrq*1e-6; % Common Gaussian dampening [Hz]
lorentzLB   = exT2 * scalingT2; % Common Lorentzian dampening [Hz]
freqShift   = 0; % Common frequency shift [Hz]
% Amplitudes for each metabolite basis function and each spline basis
% function:
ampl        = [zeros(nMets,1); zeros(nSplines,1)];

% Concatenate together into one large starting value vector.
x0          = [ph0; ph1; gaussLB; freqShift; ampl];


%%% 3. CALL NON-LINEAR SOLVER %%%
% Run the non-linear solver to optimize the non-linear parameters. At each
% iteration of the non-linear least squares solver, the linear parameters
% are calculated separately (similar to the VARiable PROjection algorithm).
%
% We're using a Levenberg-Marquardt implementation for the non-linear 
% problem that allows us to impose hard box constraints on the non-linear
% parameters, keeping them within reasonable limits.
% (c) Alexander Drentler
% https://www.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox

% Pack everything up into structs to pass it on to the solver.
% ... data:
inputData.dataToFit     = dataToFit;
inputData.resBasisSet   = reducedBasisSet;
inputData.splineArray   = splineArray;
% ... settings:
inputSettings.fitRangePPM   = fitRangePPM;
inputSettings.lorentzLB     = lorentzLB;

% Set the hard box constraints for the parameters
lb_ph0          = -15; 
ub_ph0          = +15; % Zero order phase shift [deg]
lb_ph1          = -10; 
ub_ph1          = +10;  % First order phase shift [deg/ppm]
lb_gaussLB      = 0; 
ub_gaussLB      = sqrt(5000); % Gaussian dampening [Hz]
lb_freqShift    = -15; 
ub_freqShift    = +15;    % Frequency shift [Hz]
lb_amplMets     = zeros(nMets,1);
ub_amplMets     = Inf*ones(nMets,1); % Metabolite amplitudes
lb_beta_j       = -Inf*ones(nSplines,1);
ub_beta_j       = +Inf*ones(nSplines,1); % Baseline spline amplitudes

% Concatenate together into LB/UB vectors
lb  = [lb_ph0; lb_ph1; lb_gaussLB; lb_freqShift; lb_amplMets; lb_beta_j];
ub  = [ub_ph0; ub_ph1; ub_gaussLB; ub_freqShift; ub_amplMets; ub_beta_j];

% Set up and run the non-linear solver.
%%%%%%% Old call for Dentler LM start
    % opts.Display    = 'off';
    % opts.TolFun     = 1e-6;
    % opts.TolX       = 1e-6;
    % opts.MaxIter    = 400;
    % [x] = LevenbergMarquardt(@(x) fit_Osprey_PrelimReduced_Model(x, inputData, inputSettings),x0,lb,ub,opts);
%%%%%%% Old call for Dentler LM end

inputSettings.NoiseSD = dataToFit.NoiseSD;

opts = optimoptions('lsqnonlin', ...
                    'Algorithm','levenberg-marquardt', ...      % Use LM
                    'SpecifyObjectiveGradient',true,...        % Use analytic jacobian
                    'CheckGradients',false, ...                 % Check gradient
                    'FiniteDifferenceType','central', ...       % for numerically calculated jacobian only
                    'MaxIterations',1000, ...                   % Iterations
                    'Display','none');                       % Display no iterations
[x,~,~,~,~,~,~] = lsqnonlin(@(x) fit_Osprey_PrelimReduced_Model(x, inputData, inputSettings), x0, lb, ub, opts ); % Run solver


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes and
% baseline parameters).
% [fitParamsFinal] = fit_Osprey_PrelimReduced_finalLinear(x, inputData, inputSettings);


fitParamsFinal.ph0          = x(1);
fitParamsFinal.ph1          = x(2);
fitParamsFinal.gaussLB      = x(3);
fitParamsFinal.lorentzLB    = lorentzLB;
fitParamsFinal.freqShift    = x(4);
fitParamsFinal.ampl         = x(5:4+nMets);
fitParamsFinal.beta_j       = x(5+nMets:end);


%%% 5. CREATE OUTPUT %%%
% Return the fit parameters from the final linear computation to be used in
% the next LCModel analysis step (i.e. the fit with the full basis set):
fitParamsStep1  = fitParamsFinal;


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,J] = fit_Osprey_PrelimReduced_Model(x, inputData, inputSettings)
%   This function receives the current set of non-linear parameters during
%   every iteration of the non-linear least-squares optimization. The
%   parameters are applied to the basis set.
%
%   Then, a limited-memory Broyden-Fletcher-Goldfarb-Shanno solver allowing
%   for boxed constraints is applied to determine the optimal amplitudes of 
%   the linear parameters to the current set of non-linear parameters.
%
%   The linear parameters are usually the amplitudes for the basis
%   functions and cubic baseline splines.
%
%   This function returns the difference between the input data and
%   the current optimized model. This difference is used by the NLLS
%   solver.
%
%   USAGE:
%       F = fit_LCModel_PrelimReduced_Model(x, inputData, inputSettings)
%
%   INPUTS:
%       x           = Vector providing the last set of parameters coming
%                       out of the non-linear solver
%       inputData   = Struct containing all necessary data to prepare the
%                       final linear solving step (data, basis set...).
%       inputSettings = Struct containing all necessary settings.
%
%   OUTPUTS:
%       F           = Difference between data and model. This is the
%                       objective to be least-squares-optimized by the 
%                       non-linear solver.


%%% 1. UNPACK THE INPUT %%%
% ... data:
dataToFit     = inputData.dataToFit;
resBasisSet   = inputData.resBasisSet;
splineArray   = inputData.splineArray;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
lorentzLB     = inputSettings.lorentzLB;
% ... fit parameters
nMets       = resBasisSet.nMets; % number of metabolite basis functions
nSplines    = size(splineArray,2); % number of spline basis functions
ph0   = x(1) * pi/180; % zero-order phase correction [convert from deg to rad]
ph1  = x(2) * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = x(3); % Gaussian dampening [Hz]
freqShift   = x(4); % Common frequency shift [Hz]
ampl        = x(5:end); % Amplitudes for each basis function


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
t = resBasisSet.t;
for ii=1:nMets
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) .* exp(-1i*freqShift.*t)' .* exp(-lorentzLB.*t)' .* exp(-gaussLB^2.*t.*t)';    
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) * exp(1i*ph0);
end
resBasisSet.specs = fftshift(fft(resBasisSet.fids,[],1),1);
resBasisSetWithTDOps = resBasisSet;

% Run the frequency-domain operations on the basis functions
% (first order phase correction)
% Cut out the frequency range of the basis set
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end),size(splineArray,1));
% Create a ppm vector around a pivot point (water)
ppm_ax      = resBasisSet.ppm;
pivotPoint  = 4.68;
multiplier  = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nMets
    resBasisSet.specs(:,ii) = resBasisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
resBasisSet.fids = ifft(fftshift(resBasisSet.specs,1),[],1);

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = B  * exp(1i*ph0);
B = B .* exp(1i*ph1*multiplier);
B = [real(B); imag(B)];


%%% 3. SET UP AND CALL SOLVER FOR LINEAR PARAMETERS %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Set up the equation system to be solved
% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
A   = [real(resBasisSet.specs); imag(resBasisSet.specs)];
AB  = [A B];
% Cut out the data over the fit range, and turn complex into real problem
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),size(splineArray,1));
data        = [real(dataToFit.specs); imag(dataToFit.specs)];
b           = data;

%%% 4. CREATE OBJECTIVE FUNCTION
% The objective function to be minimized by the non-linear least squares
% solver is the fit residual:
F = data - (AB*ampl);

Sigma = inputSettings.NoiseSD;                           % Get sigma
F = F./Sigma;


% Plot (comment out if not debugging)
% figure(99)
% plot(data); hold;
% plot(AB*ampl);
% plot(B*ampl(size(A,2)+1:end)); plot(data - (AB*ampl) + 1.1*max(data));
% for rr = 1:nMets
%     plot(ampl(rr)*A(:,rr));
% end
% title('Preliminary Analysis with reduced basis set');
% hold;

%%% 5. CALCULATE ANALYTIC JACOBIAN 

Acomp = resBasisSet.specs;
Bcomp = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
Bcomp = Bcomp  * exp(1i*ph0);
Bcomp = Bcomp .* exp(1i*ph1*multiplier);

ABcomp = [Acomp Bcomp];
completeFit = ABcomp*ampl;


%Computation of the Jacobian
J = zeros(length(data)/2,length(x));


%derivative wrt ph0
J(:,1) = 1i * completeFit * pi/180;

%derivative wrt ph1
J(:,2) = 1i *  completeFit .* multiplier * pi/180;

tempBasis = resBasisSetWithTDOps;
%derivative wrt gaussLB
for ii = 1:nMets
        tempBasis.fids(:,ii) = -2*gaussLB*tempBasis.fids(:,ii).*(t.^2)';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nMets
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nMets
    J(:,3) = J(:,3) + tempBasis.specs(:,ii)*ampl(ii);             
end


tempBasis = resBasisSetWithTDOps;
%derivative wrt freqShift
for ii = 1:nMets
        tempBasis.fids(:,ii) = 1i*tempBasis.fids(:,ii).*t';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nMets
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nMets
    J(:,4) = J(:,4) + tempBasis.specs(:,ii)*ampl(ii);             
end

% derivative  wrt basis set  amplitudes 
for ii=1:nMets
        J(:,4+ii) = Acomp(:,ii);
end


% derivative wrt spline amplitudes
for ii=1:nSplines
    J(:,4+nMets+ii) = Bcomp(:,ii);
end

Sigma = inputSettings.NoiseSD;                           % Get sigma
J = J./Sigma;

J   = [real(J); imag(J)]; 

J = -J;

end


