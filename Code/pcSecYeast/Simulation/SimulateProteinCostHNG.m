function SimulateProteinCostHNG(a,b)
% This function calculates the protein and ER synthesis costs under humanized glycosylation conditions

% % SCE
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecYeast/');

[~,~,protein_info] = xlsread('protein_information_SCE.xlsx');

NG = cell2mat(protein_info(2:end,6));
protein = protein_info(2:end,2);
if any(NG ~= 0)
    gene_idx = find(NG ~= 0);
    NG_protein = protein(gene_idx,1);
else
    NG_protein = [];
end

ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strrep(ERprotein,'-','');
ERprotein = strcat(ERprotein,'new');
protein_info(:,2) = strcat(protein_info(:,2),'new');
TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
ERprotein = [ERprotein;TP];
allName_NG = cell(0,1);
% load model and param
load('pcSecYeast_NG.mat','model');
load('enzymedata_NG_SCE.mat','enzymedata');
load('enzymedataSEC_NG_SCE.mat','enzymedataSEC');
load('enzymedataMachine_NG_SCE.mat','enzymedataMachine');
load('enzymedataDummyER_NG_SCE.mat','enzymedataDummyER');

%% Set model
% set medium
model = setMedia(model,2);
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production
rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
model.ub(rxn) = 0;

%% Set optimization
tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id
f = tot_protein * f_modeled_protein;
f_mito = 0.1;
f_unmodelER = tot_protein *0.046;
f_erm = 0.083;
clear tot_protein f_modeled_protein;
factor_k = 1; % doesn't matter this is to tune the saturation
mulist = [0.05 0.1];
ratios = [5E-7,1E-6,5E-6,1E-5,2E-5];

for n = a:b
    n
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
        model_tmp = changeRxnBounds(model_tmp,'r_2111',mu,'b');
        for m = 1:length(ratios)
            ratio = ratios(m);
            model_tmp = changeRxnBounds(model_tmp,[ERprotein{n},' exchange'],ratio,'b');
        
        % Set optimization
            rxnID = 'r_1714'; %minimize glucose uptake rate
            osenseStr = 'Maximize';
            name = [strrep(ERprotein{n},'new',''),'_',num2str(mu*100),'_ratio',num2str(1/ratio),'_HNG_SCE'];
            fileName = writeLPSCE(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedataNG_all,factor_k,name); 
            allName_NG = [allName_NG;{fileName}];
        end
    end
end

writeclusterfileLP(allName_NG(:),'sub_SCE')
end