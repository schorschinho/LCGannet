function [relativeCRLB,names] = fit_Osprey_CRLB(fitParams, resBasisSet,MM3coModel);
% Calculates the CRLBs from the analytical jacobian


nBasisfct = size(fitParams.ampl,1);

J = fitParams.J;

fisher = real(J'*J);                  % calculate the fisher matrix  
invFisher = pinv(fisher);                           % invert fisher matrix
crlbs = sqrt(diag(invFisher));                      % get raw CRLBs values
relativeCRLB = crlbs(4+ 2*nBasisfct : 3+ 3*nBasisfct)./fitParams.ampl(:) * 100; % Relative CRLBs for amplitudes

[relativeCRLB,names] = fit_Osprey_CombinedCRLB(invFisher, fitParams, resBasisSet.name, 4+ 2*nBasisfct,relativeCRLB,MM3coModel,resBasisSet.names);

end