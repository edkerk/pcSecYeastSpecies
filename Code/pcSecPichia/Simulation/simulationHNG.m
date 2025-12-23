%% simulationHNG

% % PP
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecPichia/');

allName_NG = cell(0,1);

load('pcSecPichia_NG.mat','model');
load('enzymedata_NG_PP.mat','enzymedata');
load('enzymedataSEC_NG_PP.mat','enzymedataSEC');
load('enzymedataMachine_NG_PP.mat','enzymedataMachine');
load('enzymedataDummyER_NG_PP.mat','enzymedataDummyER');
[~,~,protein_info] = xlsread('Protein_Information_PP.xlsx');

NG = cell2mat(protein_info(2:end,6));
protein = protein_info(2:end,2);
if any(NG ~= 0)
    gene_idx = find(NG ~= 0);
    NG_protein = protein(gene_idx,1);
else
    NG_protein = [];
end

%% Set model
% set medium
model = setMediaPP(model,2);
% set carbon source
model = changeRxnBounds(model,'Ex_glc_D',-1000,'l');% glucose
model = changeRxnBounds(model,'Ex_glyc',0,'l');
model = changeRxnBounds(model,'BIOMASS_glyc',0,'b');
model = changeRxnBounds(model,'Ex_meoh',0,'l');
model = changeRxnBounds(model,'BIOMASS_meoh',0,'b');
% set oxygen
model = changeRxnBounds(model,'Ex_o2',-1000,'l');
rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
model.ub(rxn) = 0;    
%% Set optimization
rxnID = 'Ex_glc_D'; %minimize glucose uptake rate
osenseStr = 'Maximize';
tot_protein = 0.37;
f_modeled_protein = extractModeledprotein(model,'BIOMASS','PROTEIN[c]'); %g/gProtein 
% BIOMASS is pseudo_biomass_rxn_id in the GEM
% PROTEIN[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.04;
clear tot_protein f_modeled_protein;
factor_k = 1;

%% Solve LPs
mu_list = [0.01:0.01:0.35];
ratios = [0.1:0.1:1];
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

for m = 1:length(ratios)
    ratio = ratios(m);
    enzymedataNG_all = enzymedata_all;
    
    NGgene = cell(0,1);
    for i = 1:size(enzymedata_all.subunit, 1)
        subunit_genes = enzymedata_all.subunit(i, :);
        is_gene_present = any(ismember(subunit_genes, NG_protein));
        if is_gene_present
            NGgene = [NGgene,subunit_genes];
            enzymedataNG_all.kcat(i) = ratio * enzymedata_all.kcat(i);    
        end
    end
    
    allName = cell(0,1);
    
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model_tmp = changeRxnBounds(model,'BIOMASS',mu,'b');
        disp(['mu = ' num2str(mu)]);
        name = ['NG_PP',ratio,'_',num2str(mu*100)];
        fileName  = writeLPGlc(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedataNG_all,factor_k,name);
        allName{end+1} = fileName;
    end
     
    allName_NG = [allName_NG, allName];
end

writeclusterfileLP(allName_NG(:),'sub_PP')
