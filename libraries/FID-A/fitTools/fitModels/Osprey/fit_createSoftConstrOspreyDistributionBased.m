function [penalty, penaltyJac] = fit_createSoftConstrOspreyDistributionBased(basisSet,ampl,nParams,nBasisFcts,soft,constTarget)
if nargin < 5
    constr1.met     = {'Lip13', 'Lip13', 'MM09', 'MM09', 'MM09', 'MM09', 'NAA', 'MM09', 'MM09', 'MM09', 'MM09'};
    constr1.wght    = [1        1        1       1       1       1       1          1       1       1       1];
    constr2.met     = {'Lip09', 'Lip20', 'MM20', 'MM12', 'MM14', 'MM17', 'NAAG', 'MM37', 'MM38', 'MM40', 'MM42'};
    constr2.wght    = [0.267    0.15     1.5     0.3     0.75    0.375   0.15      0.89   0.1     1.39      0.3];
    constr2.sd      = [0.1      0.07     0.375   0.1     0.45    0.1     0.15      0.40   0.1     0.375     0.1];   
else

    constr1.met     = {'Lip13', 'Lip13', 'MM09', 'MM09', 'MM09', 'MM09', 'NAA', 'MM09', 'MM09', 'MM09', 'MM09'};
    constr1.wght    = [1        1        1       1       1       1       1          1       1       1       1];
    constr2.met     = {'Lip09', 'Lip20', 'MM20', 'MM12', 'MM14', 'MM17', 'NAAG', 'MM37', 'MM38', 'MM40', 'MM42'};
    constr2.wght    = [0.267    0.15     1.5     0.3     0.75    0.375   0.15      0.89   0.1     1.39      0.3];
    constr2.sd      = [0.1      0.07     0.375   0.1     0.45    0.1     0.15      0.40   0.1     0.375     0.1];
    switch constTarget
        case 'MM'
            constr1.met{end+1} = 'MM09';               
        case 'GABA'
            constr1.met{end+1} = 'GABA'; 
    end
    constr1.wght(end+1) = 1;  
    constr2.met{end+1} = 'MM3co';
    constr2.wght(end+1) = soft;  
    constr2.sd(end+1) = 0.15; 
end



    % constr1.met     = {'NAA'};
    % constr1.wght    = [    1 ];
    % constr2.met     = {'NAAG'};
    % constr2.wght    = [0.15  ];
    % constr2.sd      = [0.15 ];

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
constr2.sd      = constr2.sd(logical(idx_toKeep));

fixIdx = zeros(size(constr1.met));
adjIdx = zeros(size(constr1.met));
fixValue = zeros(size(constr1.met));
adjValue = zeros(size(constr1.met));
for rr = 1:length(constr1.met)
    fixIdx(rr) = find(strcmpi(metsInBasisSet,constr1.met(rr)));
    adjIdx(rr) = find(strcmpi(metsInBasisSet,constr2.met(rr)));
    fixValue(rr) = ampl(fixIdx(rr));
    adjValue(rr) = ampl(adjIdx(rr));
end
actualValue = adjValue./fixValue;
actualValue(isnan(actualValue)) = 0;
actualValue(isinf(actualValue)) = 0;
expectationValue = constr2.wght;
standardDeviation = constr2.sd;

diffVec     = actualValue-expectationValue;
sdVec       = standardDeviation;


penalty  = zeros(1,nParams);
% Do correct indexing
% 3 pars for ph0,ph1,Gauss
% 2 x nBasisfunc for freqshift and lorentzian
penalty(3 +2*nBasisFcts + adjIdx) = diffVec./sdVec;

penaltyJac  = zeros(nParams,nParams);

for m_adj = 1 : length(constr1.met)
    if fixValue(m_adj) ~= 0
        penaltyJac(3 +2*nBasisFcts + adjIdx(m_adj),3 +2*nBasisFcts + fixIdx(m_adj)) = -adjValue(m_adj)./(sdVec(m_adj).*fixValue(m_adj).*fixValue(m_adj));
        penaltyJac(3 +2*nBasisFcts + adjIdx(m_adj),3 +2*nBasisFcts + adjIdx(m_adj)) = 1./(sdVec(m_adj).*fixValue(m_adj)); 
    end
end

end