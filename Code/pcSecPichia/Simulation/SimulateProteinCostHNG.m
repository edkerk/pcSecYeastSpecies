%% SimulateProteinCostHNG
% This function calculates the protein and ER synthesis costs under humanized glycosylation conditions
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecPichia/');
    
mkdir('SimulateProteinCost')
cd SimulateProteinCost

[~,~,protein_info] = xlsread('protein_information_PP.xlsx');

NG = cell2mat(protein_info(2:end,6));
protein = protein_info(2:end,2);
if any(NG ~= 0)
    gene_idx = find(NG ~= 0);
    NG_protein = protein(gene_idx,1);
else
    NG_protein = [];
end
ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strcat(ERprotein,'new');
protein_info(:,2) = strcat(protein_info(:,2),'new');
TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
TP = "PAS_chr1_" + TP;
ERprotein = [ERprotein;TP];

allName_NG = cell(0,1);
% load model and param
load('pcSecPichia_NG_PP.mat','model');
load('enzymedata_NG_PP.mat','enzymedata');
load('enzymedataSEC_NG_PP.mat','enzymedataSEC');
load('enzymedataMachine_NG_PP.mat','enzymedataMachine');
load('enzymedataDummyER_NG_PP.mat','enzymedataDummyER');

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
tot_protein = 0.37;%g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'BIOMASS','PROTEIN[c]'); %g/gProtein
% BIOMASS is pseudo_biomass_rxn_id in the GEM
% PROTEIN[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.04;
clear tot_protein f_modeled_protein;
factor_k = 1;

%% Solve LPs
mulist = [0.05 0.1];
ratios = [5E-7,1E-6,5E-6,1E-5,2E-5];

for n = a:b
    [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein(n),true); 
     % NOTE: For humanized glycosylation simulations, adjust protein MW in addTargetProtein()/calculateMW().
    [enzymedataTP] = SimulateRxnKcatCoef(model_tmp,enzymedataSEC,enzymedataTP);
    enzymedata_new = enzymedata;
    enzymedata_new.proteins = [enzymedata_new.proteins;enzymedataTP.proteins];
    enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs;enzymedataTP.proteinMWs];
    enzymedata_new.proteinLength = [enzymedata_new.proteinLength;enzymedataTP.proteinLength];
    enzymedata_new.proteinExtraMW = [enzymedata_new.proteinExtraMW;enzymedataTP.proteinExtraMW];
    enzymedata_new.kdeg = [enzymedata_new.kdeg;enzymedataTP.kdeg];
    enzymedata_new.proteinPST = [enzymedata_new.proteinPST;enzymedataTP.proteinPST];
    enzymedata_new.rxns = [enzymedata_new.rxns;enzymedataTP.rxns];
    enzymedata_new.rxnscoef = [enzymedata_new.rxnscoef;enzymedataTP.rxnscoef];
    enzymedata_new = CombineEnzymedata(enzymedata_new,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
    enzymedataNG_all = enzymedata_new;

    NGgene = cell(0,1);
    for k = 1:size(enzymedata_new.subunit, 1)
        subunit_genes = enzymedata_new.subunit(k, :);
        is_gene_present = any(ismember(subunit_genes, NG_protein));
        if is_gene_present
            NGgene = [NGgene,subunit_genes];
            enzymedataNG_all.kcat(k) = 0.8 * enzymedata_new.kcat(k);
        end
    end  
    for j = 1:length(mulist)
        disp(num2str(j))
        mu = mulist(j);
        model_tmp = changeRxnBounds(model_tmp,'BIOMASS',mu,'b');
        for m = 1:length(ratios)
            ratio = ratios(m);
            model_tmp = changeRxnBounds(model_tmp,[ERprotein{n},' exchange'],ratio,'b');
            rxnID = 'Ex_glc_D'; %minimize glucose uptake rate    
            osenseStr = 'Maximize';
            name = [strrep(ERprotein{n},'new',''),'_',num2str(mu*100),'_ratio',num2str(1/ratio),'_HNG_PP'];
            fileName = writeLPGlc(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedataNG_all,factor_k,name);% no extra constraint
            allName_NG = [allName_NG;{fileName}];
        end
    end
end

