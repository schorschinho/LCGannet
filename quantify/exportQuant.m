function exportQuant(MRSCont,saveDestination)

if ~strcmp(MRSCont.opts.fit.method, 'LCModel')
    quants = {'amplMets','tCr','rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled','AlphaCorrWaterScaled','AlphaCorrWaterScaledGroupNormed','CRLB'};
else
    quants = {'amplMets','tCr','rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled','AlphaCorrWaterScaled','AlphaCorrWaterScaledGroupNormed','CRLB','h2oarea'};
end

for ss = 1:size(MRSCont.quantify.names.metab,2)
    for q = 1 : length(quants)
        if isfield(MRSCont.quantify.tables.metab, quants(q))
            for mm = 1 : size(MRSCont.quantify.names.metab,1)
                if ~isempty(MRSCont.quantify.names.metab{mm,ss})  
                    if ~MRSCont.flags.isPRIAM
                        MRSCont.quantify.tables.metab.(quants{q}).Voxel_1{mm,ss} = PopulateJSON(MRSCont.quantify.tables.metab.(quants{q}).Voxel_1{mm,ss},quants{q});
                        osp_WriteBIDsTable(MRSCont.quantify.tables.metab.(quants{q}).Voxel_1{mm,ss}, [saveDestination  filesep MRSCont.quantify.names.SubSpectra{mm,ss} '_' quants{q} '_Voxel_1_Basis_' num2str(mm)]);
                    else                    
                        MRSCont.quantify.tables.metab.(quants{q}).Voxel_1{mm,ss} = PopulateJSON(MRSCont.quantify.tables.metab.(quants{q}).Voxel_1{mm,ss},quants{q});
                        osp_WriteBIDsTable(MRSCont.quantify.tables.metab.(quants{q}).Voxel_1{mm,ss}, [saveDestination  filesep MRSCont.quantify.names.SubSpectra{mm,ss} '_' quants{q} '_Voxel_1_Basis_' num2str(mm)]);
                        MRSCont.quantify.tables.metab.(quants{q}).Voxel_2{mm,ss} = PopulateJSON(MRSCont.quantify.tables.metab.(quants{q}).Voxel_2{mm,ss},quants{q});
                        osp_WriteBIDsTable(MRSCont.quantify.tables.metab.(quants{q}).Voxel_2{mm,ss}, [saveDestination  filesep MRSCont.quantify.names.SubSpectra{mm,ss} '_' quants{q} '_Voxel_2_Basis_' num2str(mm)]);
                    end
                end
            end
        end
    end
end


end

function[Table] = PopulateJSON(Table, Quant)

Table = addprop(Table, {'VariableLongNames'}, {'variable'}); % add long name to table properties

% Run thrrough quants to input relevant json info
switch Quant
    case'amplMets'
        L = 'Metabolite amplitude of %s';
        D = 'Raw amplitude of %s signal';
        u = 'arbitrary';
    case 'tCr'
        L = '%s to Cr ratio';
        D = 'Signal of %s represented as a ratio to the Cr signal';
        u = 'arbitrary';
    case 'rawWaterScaled'
        L = '%s to raw water ratio';
        D = 'Signal of %s represented as a ratio to the water peak';
        u = 'arbitrary';
    case 'CSFWaterScaled'
        L = '%s to raw water ratio, corrected for CSF fraction';
        D = 'Signal of %s as a ratio to raw water, but corrected for CSF fraction: WaterScaled/ (1 - fCSF)';
        u = 'arbitrary';
    case 'TissCorrWaterScaled'
        L = '%s molal concentration, using tissue corrected water';
        D = '%s scaled to tissue corrected water according to Gasparovic et al, Magn Reson Med 55:1219-26 (2006).';
        u = 'mol/Kg';
    case 'AlphaCorrWaterScaled'
        L = '%s alpha-corrected molal concentration, using tissue corrected water';
        D = '%s scaled to tissue corrected water according to Gasparovic et al, Magn Reson Med 55:1219-26 (2006). Alpha (cWM/cGM) correction applied according to: ConcIU/(fGM + alpha*fWM)';
        u = 'mol/Kg';
    case 'AlphaCorrWaterScaledGroupNormed'
        L = '%s alpha-corrected molal concentration, using tissue corrected water';
        D = '%s scaled to tissue corrected water according to Gasparovic et al, Magn Reson Med 55:1219-26 (2006). Alpha (cWM/cGM) correction applied according to: ConcIU/(fGM + alpha*fWM), with group normalization: (meanfGM + alpha * meanfWM) / ((fGM + alpha * fWM) * (meanfGM + meanfWM))';
        u = 'mol/Kg';
    case 'CRLB'
        L = '%s Cramer-Rao lower-bound';
        D = 'Lower bound of measurement error for %s';
        u = 'arbitrary';
    case 'h2oarea'
        L = 'Area ratio of %s to water peak';%CWDJ???
        D = 'Area of the water peak?';%CWDJ???
        u = 'arbitrary';
end

% input L, D, and u in to Table structure and parse this back to main function
for JJ=1:length(Table.Properties.VariableNames)
    Table.Properties.CustomProperties.VariableLongNames{JJ} = sprintf(L,Table.Properties.VariableNames{JJ});  
    Table.Properties.VariableDescriptions{JJ} = sprintf(D,Table.Properties.VariableNames{JJ});
    Table.Properties.VariableUnits{JJ} = u;
end

end