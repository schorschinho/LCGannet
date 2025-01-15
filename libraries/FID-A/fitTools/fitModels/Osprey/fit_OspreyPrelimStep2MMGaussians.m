function [fitParamsStep2] = fit_OspreyPrelimStep2MMGaussians(dataToFit, resBasisSet, fitRangePPM, fitParamsStep1)
%% [fitParamsStep2] = fit_OspreyPrelimStep2MM(dataToFit, resBasisSet, minKnotSpacingPPM, fitRangePPM, fitParamsStep1, refFWHM)
%   Performs the second step of the LCModel preliminary analysis
%   according to the LCModel algorithm. The algorithm is described in:
%       S.W. Provencher, "Estimation of metabolite concentrations from
%       localized in vivo NMR spectra", Magn Reson Med 30(6):672-679 (1993)
%
%   During the second step, the input spectrum is fit using the full basis
%   set allowing individual frequency shifts and Lorentzian dampening
%   for each metabolite basis function. Additionally, a lineshape
%   convolution is applied to account for deviations from the ideal
%   Voigtian lineshape (= Gaussian/Lorentzian convolution).

%   In addition, pre-defined macromolecule and lipid basis functions are 
%   added, as well as an unregularized baseline. MM/lipids and the baseline
%   are phased with the same phasing parameters as the metabolites, but the
%   lineshape convolution is not applied to them.
%
%   Input:
%       dataToFit       = FID-A data structure
%       basisSet        = FID-A basis set container
%       minKnotSpacing  = Scalar: minimum baseline knot spacing 
%                           (this is the DKNTMN parameter in LCModel)
%       fitRangePPM     = 2-element vector: fit range [ppm]
%                           (this is the range over which the difference 
%                           between spectrum and model is minimized)
%       fitParamsStep1  = Fit parameters from the first preliminary
%                           analysis
%       refFWHM         = Preliminary linewidth estimate (in ppm) that is
%                           used to determine the width of the lineshape
%                           function that the basis functions are
%                           subsequently convolved with
%   Output:
%       fitParamsStep2  = Fit parameters:
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

%%% 2. SET AND GET STARTING VALUES %%%
% Set the starting values for the non-linear parameters to be passed on
% to the Levenberg-Marquardt NLLS solving algorithm. Note that the
% amplitude parameters do not need to be initialized, as the
% non-negative linear Lawson-Hanson solver does not require starting
% values to be provided.
nBasisFcts  = resBasisSet.nMets + resBasisSet.nMM; % number of basis functions
ph0         = fitParamsStep1.ph0; % zero-order phase correction [deg]
ph1         = fitParamsStep1.ph1; % first-order phase correction [deg/ppm]
gaussLB     = 0; % Gaussian dampening [Hz]
lorentzLB   = ones(nBasisFcts,1); % Lorentzian dampening [Hz] for each basis function
freqShift   = ones(nBasisFcts,1); % Frequency shift [Hz] for each basis function
ampl        = zeros(nBasisFcts,1); % Amplitude parameters for basis functions and baseline

% Concatenate all initial guesses together into one large x0 vector.
x0          = [ph0; ph1; gaussLB; lorentzLB; freqShift; ampl];


%%% 3. CLL NON-LINEAR SOLVER %%%
% Run the non-linear solver to optimize the non-linear parameters.
% Pack everything up into structs to pass it on to the solver.
% ... data:
inputData.dataToFit     = dataToFit;
inputData.resBasisSet   = resBasisSet;
% ... settings:
inputSettings.fitRangePPM   = fitRangePPM;
inputSettings.regParameter  = regParameter;


% Set the hard box constraints for the parameters
nMets   = resBasisSet.nMets;
nMM     = resBasisSet.nMM;
lb_ph0              = -7.5; 
ub_ph0              = +7.5; % Zero order phase shift [deg]
lb_ph1              = -2.5; 
ub_ph1              = +2.5; % First order phase shift [deg/ppm]
lb_gaussLB          = 0; 
ub_gaussLB          = sqrt(5000); % Gaussian dampening [Hz]
lb_lorentzLB_MM     = zeros(nMM, 1);  
ub_lorentzLB_MM     =  100 * ones(nMM, 1); % Lorentzian dampening [Hz] - MM/Lipids
lb_freqShift_MM     = -6.5 * ones(nMM,1); 
ub_freqShift_MM     = +6.5 * ones(nMM,1); % Frequency shift [Hz] - MM/Lipids
lb_ampl             = zeros(nMets+nMM,1); 
ub_ampl             = +Inf * ones(nMets+nMM,1); % Amplitude for metabolite and spline basis functions

% Concatenate together into LB/UB vectors
lb = [lb_ph0; lb_ph1; lb_gaussLB; lb_lorentzLB_MM; lb_freqShift_MM; lb_ampl;];
ub = [ub_ph0; ub_ph1; ub_gaussLB; ub_lorentzLB_MM; ub_freqShift_MM; ub_ampl;];


inputSettings.NoiseSD = dataToFit.NoiseSD;

opts = optimoptions('lsqnonlin', ...
                    'Algorithm','levenberg-marquardt', ...      % Use LM
                    'SpecifyObjectiveGradient',true,...        % Use analytic jacobian
                    'CheckGradients',true, ...                 % Check gradient
                    'FiniteDifferenceType','central', ...       % for numerically calculated jacobian only
                    'MaxIterations',1000, ...                   % Iterations
                    'FunctionTolerance',1e-6,...                1e-6
                    'OptimalityTolerance',1e-4,...              1e-4
                    'StepTolerance', 1e-6,...                   1e-6
                    'Display','none');                       % Display no iterations

inputSettings.sc = 0;                                           % Don't use distribution based soft constraints 
[x,~,~,~,~,~,~] = lsqnonlin(@(x) fit_Osprey_PrelimStep2_ModelMM(x, inputData, inputSettings), x0, lb, ub, opts ); % Run solver

inputSettings.sc = 1;                                         % Run final iteration with distribution based soft constraints
lb(1:2*nBasisFcts+3) = x(1:2*nBasisFcts+3);
lb(3*nBasisFcts+4+size(splineArray,2):end) = x(3*nBasisFcts+4+size(splineArray,2):end);
x(2*nBasisFcts+4:3*nBasisFcts+3+size(splineArray,2))=zeros(nBasisFcts+nSplines,1);
ub(1:2*nBasisFcts+3) = x(1:2*nBasisFcts+3);
ub(3*nBasisFcts+4+size(splineArray,2):end) = x(3*nBasisFcts+4+size(splineArray,2):end);
[x,~,~,~,~,~,J] = lsqnonlin(@(x) fit_Osprey_PrelimStep2_Model(x, inputData, inputSettings), x, lb, ub, opts ); % Run solver with soft constraints
fitParamsFinal.ph0          = x(1);
fitParamsFinal.ph1          = x(2);
fitParamsFinal.gaussLB      = x(3);
fitParamsFinal.lorentzLB    = x(4:nBasisFcts+3);
fitParamsFinal.freqShift    = x(nBasisFcts+4:2*nBasisFcts+3);
fitParamsFinal.ampl         = x(2*nBasisFcts+4:3*nBasisFcts+3);
fitParamsFinal.beta_j       = x(3*nBasisFcts+4:3*nBasisFcts+3+size(splineArray,2));
fitParamsFinal.lineShape    = x(3*nBasisFcts+4+size(splineArray,2):end);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes and
% baseline parameters).
% [fitParamsFinal,J] = fit_Osprey_PrelimStep2_finalLinearMM(x, inputData, inputSettings);
% Return the final fit parameters

fitParamsFinal.J = J;

%%% 5. CREATE OUTPUT %%%
% Return the fit parameters from the final linear computation to be used in
% the next LCModel analysis step (i.e. the fit with the full basis set):
fitParamsStep2  = fitParamsFinal;



end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LCModel preliminary analysis step 2 model
function [F,J] = fit_Osprey_PrelimStep2_ModelMM(x, inputData, inputSettings)
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
%       F = fit_LCModel_PrelimStep2_Model(x, inputData, inputSettings)
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
stdNoise      = inputData.stdNoise;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
regParameter  = inputSettings.regParameter;

% ... fit parameters
nMets       = resBasisSet.nMets;
nMM         = resBasisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
ph0         = x(1) * pi/180; % zero-order phase correction [convert from deg to rad]
ph1         = x(2) * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = x(3); % Gaussian dampening [Hz^2]
lorentzLB   = x(4:nBasisFcts+3); % Lorentzian dampening [Hz] for each basis function
freqShift   = x(nBasisFcts+4:2*nBasisFcts+3); % Frequency shift [Hz] for each basis function
ampl        = x(2*nBasisFcts+4:3*nBasisFcts+3); % Amplitudes


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
t = resBasisSet.t;
for ii=1:nBasisFcts
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';    
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) * exp(1i*ph0);
end
resBasisSet.specs = fftshift(fft(resBasisSet.fids,[],1),1);

% Run the frequency-domain operations on the basis functions
% (first order phase correction)
% Cut out the frequency range of the basis set
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end));
% Create a ppm vector around a pivot point (water)
ppm_ax = resBasisSet.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    resBasisSet.specs(:,ii) = resBasisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
resBasisSet.fids = ifft(fftshift(resBasisSet.specs,1),[],1);

%%% 3. SET UP THE LINEAR SOLVER %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)
A = real(resBasisSet.specs);
% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
AB = [A];
% Cut out the data over the fit range, and use real part only
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end));
data        = real(dataToFit.specs);
b           = data;

if inputSettings.sc
    [penalty,penaltyJac] = fit_createSoftConstrOspreyDistributionBased(resBasisSet,ampl,length(x),nBasisFcts);
else
    penalty =[];
    penaltyJac = [];
end

%%% 4. CREATE OBJECTIVE FUNCTION
% The objective function to be minimized by the non-linear least squares
% solver is the fit residual.
%

% Return the loss function
% F = 1/sqrt(stdNoise) * SOS + regB + penaltyLBFS; 
% This would be the classic LCModel function to be minimized
% For the preliminary step, just return the functional without any regularization
F = (data - AB*ampl);

Sigma = inputSettings.NoiseSD;                           % Get sigma
F = F./Sigma;

F = vertcat(F,penalty');

if ~isempty(GAP)
     F = vertcat(F(ppm_ax<GAP(1)),F(ppm_ax>GAP(2)));   
end
%%% 5. CALCULATE ANALYTIC JACOBIAN 

completeFit = AB*ampl;

% Plot (comment out if not debugging)
% figure(99)
% plot(data); hold;
% plot(AB*ampl);
% plot(B*ampl(size(A,2)+1:end)); plot(data - (AB*ampl) + 1.1*max(data));
% for rr = 1:(nMets+nMM)
%     plot(ampl(rr)*A(:,rr));
% end
% title('Preliminary Analysis with full basis set (unregularized)');
% hold;

%Computation of the Jacobian
J = zeros(length(data),length(x));

tempBasis = resBasisSetWithTDOps;
%derivative wrt ph0
J(:,1) = 1i * completeFit * pi/180;


tempBasis = resBasisSetWithTDOps;
%derivative wrt ph1
J(:,2) = 1i * completeFit * pi/180 .* multiplier;



tempBasis = resBasisSetWithTDOps;
%derivative wrt gaussLB
for ii = 1:nBasisFcts
        tempBasis.fids(:,ii) = -tempBasis.fids(:,ii).*(t.^2)';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nBasisFcts
        J(:,3) = J(:,3) + tempBasis.specs(:,ii)*ampl(ii);             
end

tempBasis = resBasisSetWithTDOps;
%derivative wrt lorentzLB
for ii = 1:nBasisFcts
        tempBasis.fids(:,ii) = -tempBasis.fids(:,ii).*t';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nBasisFcts
        J(:,3+ii) = J(:,3+ii) + tempBasis.specs(:,ii)*ampl(ii);             
end

tempBasis = resBasisSetWithTDOps;
%derivative wrt freqShift
for ii = 1:nBasisFcts
        tempBasis.fids(:,ii) = 1i*tempBasis.fids(:,ii).*t';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nBasisFcts
        J(:,3+nBasisFcts+ii) = J(:,3+nBasisFcts+ii) + tempBasis.specs(:,ii)*ampl(ii);             
end

% derivative  wrt basis set  amplitudes 
for ii=1:nBasisFcts
        J(:,3+2*nBasisFcts+ii) = Acomp(:,ii);
end

Sigma = inputSettings.NoiseSD;                           % Get sigma
J = J./Sigma;

J = -real(J);
J = cat(1, J, penaltyJac);                                          % Append to Jacobian
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FINAL LINEAR ITERATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitParamsFinal] = fit_Osprey_PrelimStep2_finalLinearMM(x, inputData, inputSettings)
%   This function is applied after the final iteration of the non-linear
%   solver has returned the final set of non-linear parameters.
%
%   At this point, the linear solver has to be run one last time to
%   estimate the final set of linear parameters.
%
%   The function returns all model parameters.
%
%   USAGE:
%       fitParamsFinal = fit_finalLinearSolver(x, inputData, inputSettings)
%
%   INPUTS:
%       x           = Vector providing the last set of parameters coming
%                     out of the non-linear solver
%       inputData   = Struct containing all necessary data to prepare the
%                     final linear solving step (data, basis set...).
%       inputSettings = Struct containing all necessary settings.
%
%   OUTPUTS:
%       fitParamsFinal = Set of final fit parameters


%%% 1. UNPACK THE INPUT %%%
% ... data:
dataToFit     = inputData.dataToFit;
resBasisSet   = inputData.resBasisSet;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
% ... fit parameters
nMets       = resBasisSet.nMets;
nMM         = resBasisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
ph0         = x(1) * pi/180; % zero-order phase correction [convert from deg to rad]
ph1         = x(2) * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = x(3); % Gaussian dampening [Hz^2]
lorentzLB   = x(4:nBasisFcts+3); % Lorentzian dampening [Hz] for each basis function
freqShift   = x(nBasisFcts+4:2*nBasisFcts+3); % Frequency shift [Hz] for each basis function
ampl        = x(2*nBasisFcts+4:3*nBasisFcts+3);


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
t = resBasisSet.t;
for ii=1:nBasisFcts
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';    
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) * exp(1i*ph0);
end
resBasisSet.specs = fftshift(fft(resBasisSet.fids,[],1),1);
resBasisSetWithTDOps = resBasisSet;

% Run the frequency-domain operations on the basis functions
% (first order phase correction)
% Cut out the frequency range of the basis set
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end));
ppm_ax = resBasisSet.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
for ii=1:nBasisFcts
    resBasisSet.specs(:,ii) = resBasisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
resBasisSet.fids = ifft(fftshift(resBasisSet.specs,1),[],1);



%%% 3. SET UP AND CALL SOLVER FOR LINEAR PARAMETERS %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)
A = real(resBasisSet.specs);
% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
AB = [A];
% Cut out the data over the fit range, and use real part only
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end));
data        = real(dataToFit.specs);
b           = data;

% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( AB*x - b)^2;
AtA     = AB'*AB; Ab = AB'*b;
grad    = @(x) 2*( AtA*x - Ab );
hess    = @(x) 2*AtA;

% Define bounds. The lower bounds for the metabolite/MM/lipid basis
% functions are zero. All other parameters are supposed to be unbound.
l = [zeros(nMets+nMM,1);];
u = [inf*ones(nMets+nMM,1);];


% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad, hess);

opts = optimoptions("fmincon","SpecifyObjectiveGradient",true,'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,'Display','none',...
    'HessianFcn','objective',MaxFunctionEvaluations=Inf,MaxIterations=Inf,...
    StepTolerance=1e-8,FunctionTolerance=1e-12,OptimalityTolerance=1e-12);

[ampl] =...
    fmincon(fun, ampl, [],[],[],[],l,u,[],opts);


%%% 4. ADD SOFT CONSTRAINTS ON AMPLITUDES %%%
% To impose soft constraints on the amplitudes, we can augment the problem
% with additional rows in the equation system. This is done in the function
% fit_createSoftConstrOsprey.
% (see Wilson et al., MRM 2011)
A_augB  = [AB];
b_aug   = [b];

% Now, run the L-BFGS-B algorithm again with the augmented equation system
% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( A_augB*x - b_aug)^2;
AtA     = A_augB'*A_augB; A_augb = A_augB'*b_aug;
grad    = @(x) 2*( AtA*x - A_augb );
hess    = @(x) 2*AtA;

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad, hess);

opts = optimoptions("fmincon","SpecifyObjectiveGradient",true,'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,'Display','none',...
    'HessianFcn','objective',MaxFunctionEvaluations=Inf,MaxIterations=Inf,...
    StepTolerance=1e-8,FunctionTolerance=1e-12,OptimalityTolerance=1e-12);

[ampl] =...
    fmincon(fun, ampl, [],[],[],[],l,u,[],opts);


%%% 5. CREATE OUTPUT %%%
% Return the final fit parameters
fitParamsFinal.ampl         = ampl(1:size(A,2));
fitParamsFinal.ph0          = x(1);
fitParamsFinal.ph1          = x(2);
fitParamsFinal.gaussLB      = x(3);
fitParamsFinal.lorentzLB    = x(4:nBasisFcts+3);
fitParamsFinal.freqShift    = x(nBasisFcts+4:2*nBasisFcts+3);

% Plot (comment out if not debugging)
% figure(99)
% plot(data); hold;
% plot(AB*ampl);
% plot(data - (AB*ampl) + 1.1*max(data));
% for rr = 1:(nMets+nMM)
%     plot(ampl(rr)*A(:,rr));
% end
% title('Preliminary Analysis with full basis set (unregularized)');
% hold;

%Computation of the Jacobian
J = zeros(length(data),length(x));

tempBasis = resBasisSetWithTDOps;
%derivative wrt ph0
J(:,1) = 1i * completeFit * pi/180;


tempBasis = resBasisSetWithTDOps;
%derivative wrt ph1
J(:,2) = 1i * completeFit * pi/180 .* multiplier;



tempBasis = resBasisSetWithTDOps;
%derivative wrt gaussLB
for ii = 1:nBasisFcts
        tempBasis.fids(:,ii) = -tempBasis.fids(:,ii).*(t.^2)';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nBasisFcts
        J(:,3) = J(:,3) + tempBasis.specs(:,ii)*ampl(ii);             
end

tempBasis = resBasisSetWithTDOps;
%derivative wrt lorentzLB
for ii = 1:nBasisFcts
        tempBasis.fids(:,ii) = -tempBasis.fids(:,ii).*t';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nBasisFcts
        J(:,3+ii) = J(:,3+ii) + tempBasis.specs(:,ii)*ampl(ii);             
end

tempBasis = resBasisSetWithTDOps;
%derivative wrt freqShift
for ii = 1:nBasisFcts
        tempBasis.fids(:,ii) = 1i*tempBasis.fids(:,ii).*t';        
end
tempBasis = op_freqrange(tempBasis,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = tempBasis.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    tempBasis.specs(:,ii) = tempBasis.specs(:,ii) .* exp(1i*ph1*multiplier);
end

for ii = 1:nBasisFcts
        J(:,3+nBasisFcts+ii) = J(:,3+nBasisFcts+ii) + tempBasis.specs(:,ii)*ampl(ii);             
end

% derivative  wrt basis set  amplitudes 
for ii=1:nBasisFcts
        J(:,3+2*nBasisFcts+ii) = Acomp(:,ii);
end

Sigma = inputSettings.NoiseSD;                           % Get sigma
J = J./Sigma;

J = -real(J);
J = cat(1, J, penaltyJac);                                          % Append to Jacobian
end 

