function [augmA, augmb, J_augmb] = fit_createSoftConstrOsprey(basisSet, A, b, ampl)
%% [augmA, augmb] = fit_createSoftConstrOsprey(basisSet, A, b)
%   This function augments an unconstrained non-negative least-squares
%   (NNLS) problem of the form min(||Ax-b||^2) with x >= 0 to include soft 
%   constraints on the metabolite amplitude estimates x.
%
%   The soft constraint is enforced by inserting additional lines ("penalty
%   functions") into the equation system that is subsequently passed on to 
%   the NNLS solver. By varying the regularization parameter lambda 
%   (default: 0.05), the strength of the constraint can be adjusted.
%
%   USAGE:
%       [augmA, augmb] = fit_createSoftConstrOsprey(basisSet, A, b);
%
%   INPUTS:
%       basisSet    = FID-A basis set container.
%       A           = Matrix of basis functions.
%       b           = Vector of spectral data.
%       ampl        = Vector of initial amplitude estimates.
%
%   OUTPUTS:
%       augmA       = Matrix of soft constraint lines augmenting the basis
%                       function matrix A.
%       augmb       = Vector of soft constraint lines augmenting the spectral
%                       data vector b.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-06-04)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       For details of the mathematical definition of soft constraints in 
%       the context of MRS data, please refer to the TARQUIN paper:
%       Wilson et al., Magn Reson Med 65:1-12 (2011)
%       The original implementation of penalty functions is presented in:
%       Gemperline and Cash, Anal Chem 75:4236-43 (2033)
%
%   HISTORY:
%       2019-06-04: First version of the code.

% Set up vector of soft constraints. All of the below vectors require
% identical length.
if strcmp(basisSet.names{1},'A')
    constr1.met     = {'Lip13', 'Lip13', 'MM09', 'MM09', 'MM09', 'MM09', 'NAA', 'MM09', 'MM09', 'MM09', 'MM09'};
    constr1.wght    = [1        1        1       1       1       1       1          1       1       1       1];
    constr2.met     = {'Lip09', 'Lip20', 'MM20', 'MM12', 'MM14', 'MM17', 'NAAG', 'MM37', 'MM38', 'MM40', 'MM42'};
    constr2.wght    = [0.267    0.15     1.5     0.3     0.75    0.375   0.15      0.89   0.1     1.39      0.3];

%     constr1.met     = {'MM09', 'MM09', 'MM09', 'MM09', 'MM09', 'MM09', 'NAA', 'MM09', 'MM09', 'MM09', 'MM09'};
%     constr1.wght    = [1        1        1       1       1       1       1          1       1       1       1];
%     constr2.met     = {'MM12', 'MM14', 'MM17', 'MM20', 'MM22', 'MM27', 'NAAG', 'MM30', 'MM32', 'MM37', 'MM42'};
%     constr2.wght    = [1    1     1     1     1    1   0.15     1   1     1      1];

else
    constr1.met     = {'Lip13', 'Lip13', 'MM09', 'MM09', 'MM09', 'MM09', 'NAA', 'MM09', 'MM09', 'MM09', 'MM09'};
    constr1.wght    = [1        1        1       1       1       1       1          1       1       1       1];
    constr2.met     = {'Lip09', 'Lip20', 'MM20', 'MM12', 'MM14', 'MM17', 'NAAG', 'MM37', 'MM38', 'MM40', 'MM42'};
    constr2.wght    = [0.267    0.15     1.5     0.3     0.75    0.375   0.15      0.89   0.1     1.39      0.3];

%     constr1.met     = {'MM09', 'MM09', 'MM09', 'MM09', 'MM09', 'MM09', 'NAA', 'MM09', 'MM09', 'MM09', 'MM09'};
%     constr1.wght    = [1        1        1       1       1       1       1          1       1       1       1];
%     constr2.met     = {'MM12', 'MM14', 'MM17', 'MM20', 'MM22', 'MM27', 'NAAG', 'MM30', 'MM32', 'MM37', 'MM42'};
%     constr2.wght    = [1    1     1     1     1    1   0.15     1   1     1      1];
end

len1 = length(constr1.met);
len2 = length(constr1.wght);
len3 = length(constr2.met);
len4 = length(constr2.wght);
len_unique = unique([len1 len2 len3 len4]);

if length(len_unique) > 1
    error('Soft constraint vectors must have the same length. Please edit fit_createSoftConstrOsprey.m to adjust.');
end

% Check if the metabolites corresponding to the soft constraints are both
% in the basis set. Remove from the soft constraint vectors, if not.
% Check which soft-constrained metabolites are available in the basis set.
metsInBasisSet = basisSet.name;
ia1 = ismember(constr1.met, metsInBasisSet);
ia2 = ismember(constr2.met, metsInBasisSet);
idx_toKeep = zeros(length(constr1.met),1);
for rr = 1:length(idx_toKeep)
    idx_toKeep(rr) = ia1(rr) & ia2(rr);
end
constr1.met     = constr1.met(logical(idx_toKeep));
constr1.wght    = constr1.wght(logical(idx_toKeep));
constr2.met     = constr2.met(logical(idx_toKeep));
constr2.wght    = constr2.wght(logical(idx_toKeep));

% Define the strength of the constraint. Lower values (close to 0)
% correspond to a weaker constraint; higher values (1 or even multiples of
% 1) will strongly enforce the defined constraint.
lambda = 0.05;

% Calculate the augmenting lines for the b and A variables
% (Wilson et al., MRM 2011, section 2.1.6.)
augmA = zeros(length(constr1.met)+length(constr2.met), size(A,2));
augmb = zeros(length(constr1.met)+length(constr2.met), size(b,2));
J_augmb = zeros(length(constr1.met)+length(constr2.met), size(b,2));
for rr = 1:length(constr1.met)
    idx1 = find(ismember(basisSet.name,constr1.met{rr}));
    idx2 = find(ismember(basisSet.name,constr2.met{rr}));
    r_u = constr1.wght(rr);
    r_v = constr2.wght(rr);
    a_u = ampl(idx1);
    a_v = ampl(idx2);
    augmA(2*rr-1,idx1)  = lambda * norm(A(:,idx1),1);
    augmA(2*rr,idx2)    = lambda * norm(A(:,idx2),1);
    augmb(2*rr-1)       = lambda * norm(A(:,idx1),1) * ((r_u*(a_u+a_v))/(r_u+r_v));
    augmb(2*rr)         = lambda * norm(A(:,idx2),1) * ((r_v*(a_u+a_v))/(r_u+r_v));
    J_augmb(2*rr-1)       = lambda * norm(A(:,idx1),1) * (r_u*a_v)/(r_u+r_v);
    J_augmb(2*rr)         = lambda * norm(A(:,idx2),1) * (r_v*a_u)/(r_u+r_v);
end

end