function [fitParamsStep2] = fit_OspreyPrelimStep2freeGauss(dataToFit, resBasisSet, minKnotSpacingPPM, fitRangePPM, fitParamsStep1, refFWHM)
%% [fitParamsStep2] = fit_OspreyPrelimStep2freeGauss(dataToFit, resBasisSet, minKnotSpacingPPM, fitRangePPM, fitParamsStep1, refFWHM)
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
%   This specific function adds an free gaussian to estimate the co-edited
%   MM at 3 ppm of the GABA-edited difference spectrum.
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
%       2020-08-27: First version of the code.
%

%%% 1. SET UP EXPECTATION VALUES AND SPLINE BASIS FUNCTIONS %%%
% Generate the expectation values and standard deviation values for the
% delta (1/T2) Lorentzian linebroadening and the frequency shift parameter
% in LCModel.
% See LCModel manual, Section 11.14 (p. 149)
[EXT2, SDT2, SDSH]  = fit_setExpectValues(dataToFit, resBasisSet);

% Create an array of normalized cubic baseline spline basis functions.
[splineArray]       = fit_makeSplineBasis(dataToFit, fitRangePPM, minKnotSpacingPPM);
nSplines            = size(splineArray,2);


%%% 2. SET AND GET STARTING VALUES %%%
% Set the starting values for the non-linear parameters to be passed on
% to the Levenberg-Marquardt NLLS solving algorithm. Note that the
% amplitude parameters do not need to be initialized, as the
% non-negative linear Lawson-Hanson solver does not require starting
% values to be provided.
nBasisFcts  = resBasisSet.nMets + resBasisSet.nMM; % number of basis functions
ph0         = fitParamsStep1.ph0; % zero-order phase correction [deg]
ph1         = fitParamsStep1.ph1; % first-order phase correction [deg/ppm]
gaussLB     = fitParamsStep1.gaussLB; % Gaussian dampening [Hz]
lorentzLB   = EXT2'; % Lorentzian dampening [Hz] for each basis function
freqShift   = fitParamsStep1.freqShift * zeros(nBasisFcts,1); % Frequency shift [Hz] for each basis function
ampl        = zeros(nBasisFcts+nSplines,1); % Amplitude parameters for basis functions and baseline

% Create the lineshape starting values
% Set up the convolution (see LCModel manual Sec 11.12, p. 148)
RFWHM               = 2.5; % choose slightly larger than the LCM default
convolutionRange    = RFWHM * refFWHM/2;
PPMINC              = abs(dataToFit.ppm(1) - dataToFit.ppm(2));
% Calculate number of points
N_s                 = round(convolutionRange/PPMINC);
%%% OPTION A - CREATE INITIAL GUESS FROM CR SINGLET %%%
% % Apply the linebroadening parameters to the Cr singlet in order to obtain
% % an initial lineshape model
% % Find the creatine basis function
% idx_Cr          = find(strcmp(resBasisSet.name,'Cr'));
% if isempty(idx_Cr)
%     error('No basis function with nametag ''Cr'' found! Abort!');
% end
% t = resBasisSet.t;
% CrFID = resBasisSet.fids(:,idx_Cr,1) .* exp(-gaussLB^2.*t.*t)' .* exp(-lorentzLB(idx_Cr).*t)';    
% CrSPEC = fftshift(fft(CrFID,[],1),1);
% % Cut out N_s points around the maximum
% [~, maxCr_idx] = max(abs(real(CrSPEC)));
% if mod(N_s,2) == 0
%     lineShape = real(CrSPEC(maxCr_idx-N_s/2:maxCr_idx+N_s/2));
% else
%     lineShape = real(CrSPEC(maxCr_idx-floor(N_s/2):maxCr_idx+ceil(N_s/2)));
% end
%%% OPTION B - CREATE INITIAL GUESS AS DELTA FUNCTION %%%
% Initialize as delta function
if mod(N_s,2) == 0
    N_s = N_s - 1; % make odd number
end
if N_s == 0
    N_s=1;
end

if N_s < 0
    N_s=-N_s;
end

lineShape = zeros(N_s,1);
lineShape(ceil(N_s/2)) = 1;

% Normalize
lineShape = lineShape/sum(lineShape);

%Initial MM gaussian lw
GaussMMLW = 14;

% Concatenate all initial guesses together into one large x0 vector.
x0          = [ph0; ph1; gaussLB; lorentzLB; freqShift; ampl; lineShape; GaussMMLW];


%%% 3. CLL NON-LINEAR SOLVER %%%
% Run the non-linear solver to optimize the non-linear parameters. At each
% iteration of the non-linear least squares solver, the linear parameters
% are calculated separately (similar to the VARiable PROjection algorithm).
%
% We're using a Levenberg-Marquardt implementation for the non-linear 
% problem that allows us to impose hard box constraints on the non-linear
% parameters, keeping them within reasonable limits.
% (c) Alexander Drentler
% https://www.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox

% Initial regularization parameter close to zero
regParameter = 2e-2;
% Pack everything up into structs to pass it on to the solver.
% ... data:
inputData.dataToFit     = dataToFit;
inputData.resBasisSet   = resBasisSet;
inputData.splineArray   = splineArray;
% Get an estimate for the standard deviation of the noise
% This is to calculate a Q factor, but need to normalize this noise value
% here - leave alone for now
if max(dataToFit.ppm > 10)
    noiseRange = (dataToFit.ppm > 9);
else
    noiseRange = (dataToFit.ppm < -2);
end
dataNoise   = real(dataToFit.specs(noiseRange,1));
dataNoise   = detrend(dataNoise);
stdNoise    = std(dataNoise);
inputData.stdNoise      = stdNoise;
% ... settings:
inputSettings.fitRangePPM   = fitRangePPM;
inputSettings.regParameter  = regParameter;
inputSettings.EXT2          = EXT2';
inputSettings.SDT2          = SDT2';
inputSettings.SDSH          = SDSH';

% Set the hard box constraints for the parameters
nMets   = resBasisSet.nMets;
nMM     = resBasisSet.nMM;
lb_ph0              = -45; 
ub_ph0              = +45; % Zero order phase shift [deg]
lb_ph1              = -10; 
ub_ph1              = +10; % First order phase shift [deg/ppm]
lb_gaussLB          = 0; 
ub_gaussLB          = sqrt(5000); % Gaussian dampening [Hz]
lb_lorentzLB_mets   = zeros(nMets, 1); 
ub_lorentzLB_mets   = 10.0 * ones(nMets, 1); % Lorentzian dampening [Hz] - Metabolites
lb_lorentzLB_MM     = zeros(nMM, 1); 
ub_lorentzLB_MM     =  100 * ones(nMM, 1); % Lorentzian dampening [Hz] - MM/Lipids
lb_freqShift_mets   = -4.0 * ones(nMets,1); 
ub_freqShift_mets   = +4.0 * ones(nMets,1); % Frequency shift [Hz] - Metabolites
lb_freqShift_MM     = -6.5 * ones(nMM,1); 
ub_freqShift_MM     = +6.5 * ones(nMM,1); % Frequency shift [Hz] - MM/Lipids
lb_ampl             = -Inf * ones(nMets+nMM+size(splineArray,2),1); 
ub_ampl             = +Inf * ones(nMets+nMM+size(splineArray,2),1); % Amplitude for metabolite and spline basis functions
lb_lineShape        = -Inf * ones(size(lineShape));
ub_lineShape        = +Inf * ones(size(lineShape)); % Lineshape coefficients
lb_GaussMMLW        = 1;
ub_GaussMMLW        = 20; % Gaussian co-edited MM linewidth

% Concatenate together into LB/UB vectors
lb = [lb_ph0; lb_ph1; lb_gaussLB; lb_lorentzLB_mets; lb_lorentzLB_MM; lb_freqShift_mets; lb_freqShift_MM; lb_ampl; lb_lineShape; lb_GaussMMLW];
ub = [ub_ph0; ub_ph1; ub_gaussLB; ub_lorentzLB_mets; ub_lorentzLB_MM; ub_freqShift_mets; ub_freqShift_MM; ub_ampl; ub_lineShape; ub_GaussMMLW];


% Set up and run the non-linear solver.
%%%%%%% Old call for Dentler LM start
    % opts.Display    = 'off';
    % opts.TolFun     = 1e-6;
    % opts.TolX       = 1e-6;
    % opts.MaxIter    = 400;
    % [x,Res,~,~,~,~,~,J] = LevenbergMarquardt(@(x) fit_Osprey_PrelimStep2_Model(x, inputData, inputSettings), x0, lb, ub, opts);
%%%%%%% Old call for Dentler LM end

opts = optimoptions('lsqnonlin', ...
                    'Algorithm','levenberg-marquardt', ...      % Use LM
                    'SpecifyObjectiveGradient',false,...        % Use analytic jacobian
                    'CheckGradients',false, ...                 % Check gradient
                    'FiniteDifferenceType','central', ...       % for numerically calculated jacobian only
                    'MaxIterations',1000, ...                   % Iterations
                    'Display','none');                       % Display no iterations
[x,~,Res,~,~,~,J] = lsqnonlin(@(x) fit_Osprey_PrelimStep2_Model(x, inputData, inputSettings), x0, lb, ub, opts ); % Run solver


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes and
% baseline parameters).
[fitParamsFinal,resBasisSetUpdated] = fit_Osprey_PrelimStep2_finalLinear(x, inputData, inputSettings);
fitParamsFinal.Res = Res;
fitParamsFinal.J = J;
fitParamsFinal.resBasisSetUpdated = resBasisSetUpdated;

%%% 5. CREATE OUTPUT %%%
% Return the fit parameters from the final linear computation to be used in
% the next LCModel analysis step (i.e. the fit with the full basis set):
fitParamsStep2  = fitParamsFinal;


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LCModel preliminary analysis step 2 model
function [F,J] = fit_Osprey_PrelimStep2_Model(x, inputData, inputSettings)
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
splineArray   = inputData.splineArray;
stdNoise      = inputData.stdNoise;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
regParameter  = inputSettings.regParameter;
EXT2          = inputSettings.EXT2;
SDT2          = inputSettings.SDT2;
SDSH          = inputSettings.SDSH;
% ... fit parameters
nMets       = resBasisSet.nMets;
nMM         = resBasisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
nSplines    = size(splineArray,2); % number of spline basis functions
ph0         = x(1) * pi/180; % zero-order phase correction [convert from deg to rad]
ph1         = x(2) * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = x(3); % Gaussian dampening [Hz^2]
lorentzLB   = x(4:nBasisFcts+3); % Lorentzian dampening [Hz] for each basis function
freqShift   = x(nBasisFcts+4:2*nBasisFcts+3); % Frequency shift [Hz] for each basis function
ampl        = x(2*nBasisFcts+4:3*nBasisFcts+3+nSplines); % Amplitudes
lineShape   = x(3*nBasisFcts+4+nSplines:end-1); % Lineshape coefficients
% Normalize the lineshape
lineShape = lineShape/sum(lineShape);
GaussMMLW = x(end); % co-edited MM Gaussian linewidth

%%% 2. Create co-edited MM Gaussian %%%
MM3co  = op_gaussianPeak(resBasisSet.coMM3.n,resBasisSet.coMM3.sw,resBasisSet.coMM3.Bo,resBasisSet.coMM3.centerFreq,GaussMMLW,3,resBasisSet.coMM3.MMArea);
MM3co  = op_dccorr(MM3co,'p');
MM3co.nMets = 0;
MM3co.nMM = 1;
[MM3co] = fit_resampleBasis(dataToFit, MM3co);
resBasisSet.fids(:,end) = MM3co.fids(:);
resBasisSet.specs(:,end) = MM3co.specs(:);

%%% 3. APPLY THE NON-LINEAR PARAMETERS %%%
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
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = resBasisSet.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    resBasisSet.specs(:,ii) = resBasisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
resBasisSet.fids = ifft(fftshift(resBasisSet.specs,1),[],1);

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = B  * exp(1i*ph0);
B = B .* exp(1i*ph1*multiplier);
B = [real(B)];


%%% 4. SET UP THE LINEAR SOLVER %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
A = real(resBasisSet.specs);
for kk = 1:resBasisSet.nMets
    A(:,kk) = conv(A(:,kk), lineShape, 'same');
end

% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
AB = [A B];
% Cut out the data over the fit range, and use real part only
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),length(splineArray(:,1,1)));
data        = real(dataToFit.specs);
b           = data;

% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( AB*x - b)^2;
AtA     = AB'*AB; Ab = AB'*b;
grad    = @(x) 2*( AtA*x - Ab );
hess    = @(x) 2*AtA;

% Define bounds. The lower bounds for the metabolite/MM/lipid basis
% functions are zero. All other parameters are supposed to be unbound.
l = [zeros(nMets+nMM,1);    -inf*ones(nSplines,1)];
u = [inf*ones(nMets+nMM,1);  inf*ones(nSplines,1)];

%%%%%%% Old call for lbfgsb start
% Prepare the function wrapper
    % fun     = @(x)fminunc_wrapper( x, fcn, grad);
    % % Request very high accuracy:
    % opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10, 'max_its', 750);
    % opts.printEvery     = 0;
    % 
    % % Run the algorithm:
    % % Feed initial guess from the input parameters
    % opts.x0 = ampl;
    % [ampl, ~, ~] = lbfgsb(fun, l, u, opts );
%%%%%%% Old call for lbfgsb end

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad, hess);

opts = optimoptions("fmincon","SpecifyObjectiveGradient",true,'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,'Display','none',...
    'HessianFcn','objective',MaxFunctionEvaluations=Inf,MaxIterations=750,...
    StepTolerance=1e-8,FunctionTolerance=1e-12,OptimalityTolerance=1e-12);

[ampl] =...
    fmincon(fun, ampl, [],[],[],[],l,u,[],opts);

%%% 5. ADD SOFT CONSTRAINTS ON AMPLITUDES %%%
% To impose soft constraints on the amplitudes, we can augment the problem
% with additional rows in the equation system. This is done in the function
% fit_createSoftConstrOsprey.
% (see Wilson et al., MRM 2011)
[augmA, augmb] = fit_createSoftConstrOsprey(resBasisSet, AB, b, ampl);
A_augB  = [AB; augmA];
b_aug   = [b; augmb];

% Now, run the L-BFGS-B algorithm again with the augmented equation system
% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( A_augB*x - b_aug)^2;
AtA     = A_augB'*A_augB; A_augb = A_augB'*b_aug;
grad    = @(x) 2*( AtA*x - A_augb );
hess    = @(x) 2*AtA;

%%%%%%% Old call for lbfgsb start
% Prepare the function wrapper
    % fun     = @(x)fminunc_wrapper( x, fcn, grad);
    % % Request very high accuracy:
    % opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10, 'max_its', 750);
    % opts.printEvery     = 0;
    % 
    % % Run the algorithm:
    % % Feed initial guess from the input parameters
    % opts.x0 = ampl;
    % [ampl, ~, ~] = lbfgsb(fun, l, u, opts );
%%%%%%% Old call for lbfgsb end

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad, hess);

opts = optimoptions("fmincon","SpecifyObjectiveGradient",true,'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,'Display','none',...
    'HessianFcn','objective',MaxFunctionEvaluations=Inf,MaxIterations=750,...
    StepTolerance=1e-8,FunctionTolerance=1e-12,OptimalityTolerance=1e-12);

[ampl] =...
    fmincon(fun, ampl, [],[],[],[],l,u,[],opts);


%%% 6. CREATE OBJECTIVE FUNCTION
% The objective function to be minimized by the non-linear least squares
% solver is the fit residual.
%
% The LCModel algorithm wants to minimize not only the sum of squares
% between the data and the model, but imposes regularization on several
% parameters.
% Therefore, we have to reconstruct the fit spectrum from the parameters we
% derived, and then calculate the sum of squares ourselves, in addition to
% the penalty terms.
% For the last step of the preliminary analysis, this is done without any
% regularization at all.

% Calculate the sum of squares
SOS = sum((data - AB*ampl).^2);

% Calculate the baseline regularization penalty term
% Extract the baseline coefficients
beta_j    = ampl(size(A,2)+1:end);
% Create the analytical regularizor matrix specified in Eq (3.13) in the
% original publication of the CONTIN algorithm behind LCModel:
% Provencher, Comp Phys Comm 27:213-227 (1982)
% First, create a banded Toeplitz matrix
n = length(beta_j);
e = ones(n,1);
K = spdiags([1*e 0*e -9*e 16*e -9*e 0*e 1*e], -3:3, n, n);
K = full(K);
% Then, correct the first two and last two lines according to the above
% paper
K(1, 1:4)               = [8 -6 0 1];
K(2, 1:5)               = [-6 14 -6 0 1];
K(end-1, end-4:end)     = [1 0 -6 14 -6];
K(end, end-3:end)       = [1 0 -6 8];
% Bring in the factor of 1/6
K = 1/6 .* K;
% Additionally, bring in the geometric factor delta, which is the knot
% spacing, ie delta = range/(number of segments).
% !!! THIS ISN'T CLEARLY UNDERSTOOD YET !!!
normalizationFactor = 1/(size(splineArray,1)/(size(splineArray,2)-1))^3;
% Final regularization penalty
regB = norm(regParameter*sqrtm(full(K))*beta_j*1/(mean(beta_j)))^2;

% Calculate the penalty term for the Lorentzian linebroadening and
% frequency shifts
penaltyLBFS = sum((lorentzLB - EXT2).^2./(SDT2).^2 + (freqShift.^2)./(SDSH.^2));

% Return the loss function
% F = 1/sqrt(stdNoise) * SOS + regB + penaltyLBFS; 
% This would be the classic LCModel function to be minimized
% For the preliminary step, just return the functional without any regularization
F = (data - AB*ampl);

%%% 7. CALCULATE ANALYTIC JACOBIAN 
% j = sqrt(-1); % i
% 
% % Model function
% % Y(f) = exp(i(ph0+ph1(f))) * (sum(betaj * B) + sum(amp conv(M, lineshape))
% % with M = fft(basisSet.fids .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)')
% % The derivatives are 'easier' in the time-domain, so we will stick with
% % the FIDS for the calculations
% 
% % dimensions and number of matrix entries
% [npoints,~] = size(resBasisSet.fids); %number of points and number of basis functions
% t = resBasisSet.t;
% 
% % % Apply phasing to the spline basis functions
% Acomp = resBasisSet.specs;
% for kk = 1:resBasisSet.nMets
%     Acomp(:,kk) = conv(Acomp(:,kk), lineShape, 'same');
% end
% 
% Bcomp = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
% Bcomp = Bcomp  * exp(1i*ph0);
% Bcomp = Bcomp .* exp(1i*ph1*multiplier);
% 
% ABcomp = [Acomp Bcomp];
% completeFit = ABcomp*ampl;
% 
% 
% %Computation of the Jacobian
% % The size is MxN with M (number of estimated parameters) and N (number of
% % points).
% %Store the indices of the different partial derivatives 
% ph0Col = 1;
% ph1Col = 2;
% gaussLBCol = 3 * ones(1,nMets + nMM);
% lorentzLBCol = gaussLBCol(end) + 1 : (gaussLBCol(1) + nMets + nMM);
% freqShiftCol = lorentzLBCol(end) + 1 : (lorentzLBCol(end) + nMets + nMM);
% AmplCol = freqShiftCol(end) + 1 : (freqShiftCol(end) + nMets + nMM);
% SplineAmplCol = AmplCol(end) + 1 : (AmplCol(end) + size(splineArray,2));
% lineshapeCol = SplineAmplCol(end)+1 : SplineAmplCol(end) + length(lineShape);
% 
% nparams = 3 + length(lorentzLBCol) + length(freqShiftCol) + length(AmplCol) + length(SplineAmplCol) + length(lineshapeCol);
% 
% J = zeros(npoints,nparams);
% Jfd = zeros(npoints,nparams);
% 
% 
% %derivative wrt ph0
% % for ii = 1:nMets
% %     col = ph0Col(ii);
% %     J(:,col) = J(:,col)+j*basisSet.fids(:,ii)*ampl(ii);
% %     Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
% %     Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same') + B * beta_j;    
% % end
% Jfd(:,1) = j * completeFit;
% 
% %derivative wrt ph1
% % for ii = 1:nBasisFcts
% %     col = ph1Col(ii);
% %     J(:,col) = J(:,col)+j*basisSet.fids(:,ii)*ampl(ii);
% %     Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
% %     Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same') + B * beta_j;    
% % end
% Jfd(:,2) = j *  completeFit .* multiplier;
% 
% %derivative wrt gaussLB
% for ii = 1:nBasisFcts
%     if ii <= nMets % Sum up derivarives of all metabolite functions first
%         col = gaussLBCol(ii);
%         J(:,col) = J(:,col)-resBasisSet.fids(:,ii).*(t.^2)';
%         Jfd(:,col) = Jfd(:,col) +  conv(fftshift(fft(-resBasisSet.fids(:,ii).*(t.^2)',[],1),1), lineShape, 'same')*ampl(ii); 
%     else  % No convolution is applied to the MM functions
%         J(:,col) = J(:,col)-resBasisSet.fids(:,ii).*(t.^2)';
%         Jfd(:,col) = Jfd(:,col) +  fftshift(fft(-resBasisSet.fids(:,ii).*(t.^2)',[],1),1)*ampl(ii);             
%     end
% end
% 
% 
% %derivative wrt lorentzLB
% for ii = 1:nBasisFcts
%     if ii <= nMets % Sum up derivarives of all metabolite functions first
%         col = lorentzLBCol(ii);
%         J(:,col) = J(:,col)-resBasisSet.fids(:,ii).*t';
%         Jfd(:,col) = Jfd(:,col) +  conv(fftshift(fft(J(:,col),[],1),1), lineShape, 'same')*ampl(ii);  
%     else % No convolution is applied to the MM functions
%         col = lorentzLBCol(ii);
%         J(:,col) = J(:,col)-resBasisSet.fids(:,ii).*t';
%         Jfd(:,col) = Jfd(:,col) +  fftshift(fft(J(:,col),[],1),1)*ampl(ii);         
%     end
% end
% 
% %derivative wrt freqShift
% for ii=1:nBasisFcts
%     if ii <= nMets  % Sum up derivarives of all metabolite functions first
%         col = freqShiftCol(ii);
%         J(:,col) = J(:,col)-j*resBasisSet.fids(:,ii).*t';
%         Jfd(:,col) = Jfd(:,col) +  conv(fftshift(fft(J(:,col),[],1),1), lineShape, 'same')*ampl(ii);  
%     else % No convolution is applied to the MM functions
%         col = freqShiftCol(ii);
%         J(:,col) = J(:,col)-j*resBasisSet.fids(:,ii).*t';
%         Jfd(:,col) = Jfd(:,col) + fftshift(fft(J(:,col),[],1),1)*ampl(ii);          
%     end
% end
% 
% % derivative  wrt basis set  amplitudes 
% for ii=1:nBasisFcts
%     if ii <= nMets
%         col = AmplCol(ii);
%         J(:,col) = resBasisSet.fids(:,ii);
%         Jfd(:,col) = conv(fftshift(fft(J(:,col),[],1),1), lineShape, 'same');
%     else
%         col = AmplCol(ii);
%         J(:,col) = resBasisSet.fids(:,ii);
%         Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);    
%     end
% end
% 
% 
% % derivative wrt spline amplitudes
% for ii=1:length(SplineAmplCol)
%     col = SplineAmplCol(ii);
%     Jfd(:,col) = Jfd(:,col)+Bcomp(:,ii);
% end
% 
% 
% %derivative wrt lineshape 
% % We will do a discrete convolution of S'*M by using the Toeplitz matrix
% % form of the partial derivative of the lineshape vector S and M beeing 
% % the metabolite basis functions. This is allowed as convolutions are 
% % commutative and (M*S)' = M'*S = M*S' -> (M*S)' = S'*M.
% 
% %We need to create S' as a Toeplitz matrix. The derivatives will
% %essantially be ones on the diagonal (first lineshape coeff) or the upper off
% %diagonal (all other lineshape coeff). 
% 
% nLineShape = length(lineShape);
% nPoints = length(resBasisSet.fids(:,1));
% 
% %Set up the first rows of each partial derivative
% Toep1row = zeros(nLineShape,nPoints);
% for ii = 1 : nLineShape
%     Toep1row(ii,ii) = 1; 
% end
% 
% %We will set it up as a 3D vector with the partial derivatives in the third dimensions and the
% %Toeplitz matrix spanning the first two dimensions.
% ToepLineShape = zeros(nPoints,nPoints,nLineShape);
% for ii = 1 : nLineShape
%     if ii <= 1
%         ToepLineShape(:,:,ii) = toeplitz(Toep1row(ii,:));
%     else
%         ToepLineShape(:,:,ii) = tril(toeplitz(Toep1row(ii,:))); % extract lower triangle of the Toeplitz matrix
%     end
% end
% 
% for ii=1:nLineShape % Loop over the lineshape derivatives
%     col = lineshapeCol(ii);
%     for kk = 1 : nMets % Convolute all basis functions to create the full model function
%         J(:,col) = J(:,col)+resBasisSet.fids(:,kk);
%         Jfd(:,col) = Jfd(:,col) + ampl(kk)*(ToepLineShape(:,:,ii) * fftshift(fft(resBasisSet.fids(:,kk),[],1),1));
%     end
%         
% end
% 
% J = real(Jfd);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FINAL LINEAR ITERATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitParamsFinal,resBasisSetUpdated] = fit_Osprey_PrelimStep2_finalLinear(x, inputData, inputSettings)
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
splineArray   = inputData.splineArray;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
% ... fit parameters
nMets       = resBasisSet.nMets;
nMM         = resBasisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
nSplines    = size(splineArray,2);
ph0         = x(1) * pi/180; % zero-order phase correction [convert from deg to rad]
ph1         = x(2) * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = x(3); % Gaussian dampening [Hz^2]
lorentzLB   = x(4:nBasisFcts+3); % Lorentzian dampening [Hz] for each basis function
freqShift   = x(nBasisFcts+4:2*nBasisFcts+3); % Frequency shift [Hz] for each basis function
ampl        = x(2*nBasisFcts+4:3*nBasisFcts+3+size(splineArray,2));
lineShape   = x(3*nBasisFcts+4+size(splineArray,2):end-1);
% Normalize
lineShape   = lineShape/sum(lineShape);
GaussMMLW   = x(end); 

%%% 2. Create co-edited MM Gaussian %%%
MM3co  = op_gaussianPeak(resBasisSet.coMM3.n,resBasisSet.coMM3.sw,resBasisSet.coMM3.Bo,resBasisSet.coMM3.centerFreq,GaussMMLW,3,resBasisSet.coMM3.MMArea);
MM3co  = op_dccorr(MM3co,'p');
MM3co.nMets = 0;
MM3co.nMM = 1;
[MM3co] = fit_resampleBasis(dataToFit, MM3co);
resBasisSet.fids(:,end) = MM3co.fids(:);
resBasisSet.specs(:,end) = MM3co.specs(:);
resBasisSetUpdated = resBasisSet;


%%% 3. APPLY THE NON-LINEAR PARAMETERS %%%
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
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
ppm_ax = resBasisSet.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
for ii=1:nBasisFcts
    resBasisSet.specs(:,ii) = resBasisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
resBasisSet.fids = ifft(fftshift(resBasisSet.specs,1),[],1);

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = B  * exp(1i*ph0);
B = B .* exp(1i*ph1*multiplier);
B = [real(B)];


%%% 4. SET UP AND CALL SOLVER FOR LINEAR PARAMETERS %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Set up the equation system to be solved
% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
A = real(resBasisSet.specs);
for kk = 1:resBasisSet.nMets
    A(:,kk) = conv(A(:,kk), lineShape, 'same');
end

% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
AB = [A B];
% Cut out the data over the fit range, and use real part only
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),length(splineArray(:,1,1)));
data        = real(dataToFit.specs);
b           = data;

% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( AB*x - b)^2;
AtA     = AB'*AB; Ab = AB'*b;
grad    = @(x) 2*( AtA*x - Ab );
hess    = @(x) 2*AtA;

% Define bounds. The lower bounds for the metabolite/MM/lipid basis
% functions are zero. All other parameters are supposed to be unbound.
l = [zeros(nMets+nMM,1);    -inf*ones(nSplines,1)];
u = [inf*ones(nMets+nMM,1);  inf*ones(nSplines,1)];

%%%%%%% Old call for lbfgsb start
% Prepare the function wrapper
    % fun     = @(x)fminunc_wrapper( x, fcn, grad);
    % % Request very high accuracy:
    % opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10, 'max_its', 750);
    % opts.printEvery     = 0;
    % 
    % % Run the algorithm:
    % % Feed initial guess from the input parameters
    % opts.x0 = ampl;
    % [ampl, ~, ~] = lbfgsb(fun, l, u, opts );
%%%%%%% Old call for lbfgsb end

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad, hess);

opts = optimoptions("fmincon","SpecifyObjectiveGradient",true,'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,'Display','none',...
    'HessianFcn','objective',MaxFunctionEvaluations=Inf,MaxIterations=750,...
    StepTolerance=1e-8,FunctionTolerance=1e-12,OptimalityTolerance=1e-12);

[ampl] =...
    fmincon(fun, ampl, [],[],[],[],l,u,[],opts);


%%% 5. ADD SOFT CONSTRAINTS ON AMPLITUDES %%%
% To impose soft constraints on the amplitudes, we can augment the problem
% with additional rows in the equation system. This is done in the function
% fit_createSoftConstrOsprey.
% (see Wilson et al., MRM 2011)
[augmA, augmb] = fit_createSoftConstrOsprey(resBasisSet, AB, b, ampl);
A_augB  = [AB; augmA];
b_aug   = [b; augmb];

% Now, run the L-BFGS-B algorithm again with the augmented equation system
% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( A_augB*x - b_aug)^2;
AtA     = A_augB'*A_augB; A_augb = A_augB'*b_aug;
grad    = @(x) 2*( AtA*x - A_augb );
hess    = @(x) 2*AtA;

%%%%%%% Old call for lbfgsb start
% Prepare the function wrapper
    % fun     = @(x)fminunc_wrapper( x, fcn, grad);
    % % Request very high accuracy:
    % opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10, 'max_its', 750);
    % opts.printEvery     = 0;
    % 
    % % Run the algorithm:
    % % Feed initial guess from the input parameters
    % opts.x0 = ampl;
    % [ampl, ~, ~] = lbfgsb(fun, l, u, opts );
%%%%%%% Old call for lbfgsb end

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad, hess);

opts = optimoptions("fmincon","SpecifyObjectiveGradient",true,'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,'Display','none',...
    'HessianFcn','objective',MaxFunctionEvaluations=Inf,MaxIterations=750,...
    StepTolerance=1e-8,FunctionTolerance=1e-12,OptimalityTolerance=1e-12);

[ampl] =...
    fmincon(fun, ampl, [],[],[],[],l,u,[],opts);


%%% 6. CREATE OUTPUT %%%
% Return the final fit parameters
fitParamsFinal.ampl         = ampl(1:size(A,2));
fitParamsFinal.ph0          = x(1);
fitParamsFinal.ph1          = x(2);
fitParamsFinal.gaussLB      = x(3);
fitParamsFinal.lorentzLB    = x(4:nBasisFcts+3);
fitParamsFinal.freqShift    = x(nBasisFcts+4:2*nBasisFcts+3);
fitParamsFinal.lineShape    = x(3*nBasisFcts+4+size(splineArray,2):end-1);
fitParamsFinal.beta_j       = ampl(size(A,2)+1:end);
fitParamsFinal.GaussMMLW       = GaussMMLW;

% % Plot (comment out if not debugging)
% figure(99)
% plot(data); hold;
% plot(AB*ampl);
% plot(B*ampl(size(A,2)+1:end)); plot(data - (AB*ampl) + 1.1*max(data));
% for rr = 1:(nMets+nMM)
%     plot(ampl(rr)*A(:,rr));
% end
% title('Preliminary Analysis with full basis set (unregularized)');
% hold;


end 

