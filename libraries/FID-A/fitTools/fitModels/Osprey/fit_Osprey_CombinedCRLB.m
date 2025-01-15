function [crlbs,names] = fit_Osprey_CombinedCRLB(invFisher, fitParams, metaboliteNames, StartIdxMetab,crlbs)
%%  [crlbs,names] = fit_Osprey_CombinedCRLB(jacobian, metaboliteNames, parametrizations)
%   This method calculates the CRLBs of metabolite combinations inthe OspreyFitObj 
%
%   USAGE:
%       [obj] = calculateCombinedCRLB(obj,jac, metaboliteNames, parametrizations)
%
%   INPUTS:
%       invFisher   = inverse fisher matrix
%       xk          = final parameter vector
%       metaboliteNames  = metabolites in model
%       parametrizations  = struct with parametrizations in xk vector
%       
%   OUTPUTS:
%       obj       = OspreyFitObj with updated parametrizations.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)

%% 0. Set names cell arrays and houskeeping
typicalMetaboliteCombinations = {'NAA','NAAG';'GPC','PCh';'Cr','PCr';'Glu','Gln';'EA','PE';'GABA','MM3co'; 'GABA', 'MM3to2'}; 
MetaboliteCombinationNames = {'tNAA','tCho','tCr','Glx','tEA','GABAplus','GABAplus'};


nPars = 1;
pos = 1;
names = [];

BMAT = zeros(size(invFisher,1),length(metaboliteNames));                          % Some LCModel nostalgia 
ll = 1;
for kk = StartIdxMetab:StartIdxMetab + length(metaboliteNames) - 1
   BMAT(kk,ll)=1;
   ll = ll + 1;
end

%% 1. Get metabolite index, update Jacobian, update amplitude estimates
AddedMetaboliteCombinations = 0;
for mm = 1 : length(MetaboliteCombinationNames)
    idx_1 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,1}));        
    idx_2 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,2}));
    if  ~isempty(idx_1) && ~isempty(idx_2)
        AddedMetaboliteCombinations = AddedMetaboliteCombinations + 1;
        BMAT(StartIdxMetab + idx_1 - 1,end+1)=1;
        BMAT(StartIdxMetab + idx_2 - 1,end)=1;  
    end
end



%% 2. Calculate CRLB
if size(BMAT,2) > length(metaboliteNames)*nPars         % added new combinations 
    combinations = fitParams.ampl' * BMAT(StartIdxMetab:StartIdxMetab + length(metaboliteNames) - 1,:);  % build combinations of amplitude parameters
    DAPOSI = BMAT' * invFisher *BMAT;                   % multiply with inverse fisher matrix
    crlbs = sqrt(diag(DAPOSI));                         % get raw CRLBs values                    
    crlbs= (crlbs ./ combinations') * 100; % Relative CRLBs for combined amplitudes
    names = MetaboliteCombinationNames(1:AddedMetaboliteCombinations);
end
end
