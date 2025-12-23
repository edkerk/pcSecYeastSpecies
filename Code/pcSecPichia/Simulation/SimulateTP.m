%% SimulateTP

%% PP
% max production rate under various growth rates on minimal media
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecPichia/');

TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};

mulist = [0.02:0.02:0.08,0.08:0.01:0.16,0.16:0.02:0.34];

mkdir('SimulateTP_PP');
cd SimulateTP_PP/;
allname = cell(0,1);
TP = "PAS_chrx_" + TP;

%% Solve LPs
for i = 1:length(TP)
    % load model and param
    load('enzymedata_PP.mat')
    load('enzymedataMachine_PP.mat')
    load('enzymedataSEC_PP.mat')
    load('enzymedataDummyER_PP.mat')
    load('pcSecPichia.mat')


    %% Set model
    % set medium
    model = setMediaPP(model,4);
    % set carbon source
    model = changeRxnBounds(model,'Ex_glc_D',-1000,'l');
    model = changeRxnBounds(model,'BIOMASS',1000,'u');
    model = changeRxnBounds(model,'LIPIDS',1000,'u');
    model = changeRxnBounds(model,'PROTEINS',1000,'u');
    model = changeRxnBounds(model,'STEROLS',1000,'u');

    model = changeRxnBounds(model,'Ex_glyc',0,'l');
    model = changeRxnBounds(model,'BIOMASS_glyc',0,'b');
    model = changeRxnBounds(model,'LIPIDS_glyc',0,'b');
    model = changeRxnBounds(model,'PROTEINS_glyc',0,'b');
    model = changeRxnBounds(model,'STEROLS_glyc',0,'b');

    model = changeRxnBounds(model,'Ex_meoh',0,'l');
    model = changeRxnBounds(model,'BIOMASS_meoh',0,'b');
    model = changeRxnBounds(model,'LIPIDS_meoh',0,'b');
    model = changeRxnBounds(model,'PROTEINS_meoh',0,'b');
    model = changeRxnBounds(model,'STEROLS_meoh',0,'b');

    % set oxygen
    model = changeRxnBounds(model,'Ex_o2',-1000,'l');

    rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
    model.ub(rxn) = 0;
    %% Set optimization
    
    tot_protein = 0.37; % g/gCDW, total protein content (0.37 for glucose condition, 0.40 for methanol condition)
% %     tot_protein = 0.4; 

    f_modeled_protein = extractModeledprotein(model,'BIOMASS','PROTEIN[c]'); %g/gProtein
    % BIOMASS is pseudo_biomass_rxn_id in the GEM
    % PROTEIN[c] is protein id
    
    f = tot_protein * f_modeled_protein;
    f_mito = 0.1;
    f_unmodelER = tot_protein * 0.04;
    f_erm = 0.083;
    clear tot_protein f_modeled_protein;
    
    factor_k = 1; % doesn't matter this is to tune the saturation
    rxnID = strcat(TP{i},' exchange');
% %     rxnID = 'dilute_dummy';
    osenseStr = 'Maximize';
    [model_tmp,enzymedataTP] = addTargetProtein(model,TP(i));
    model = model_tmp;
    save(['modelPP_',TP{i},'.mat'],'model')
    [enzymedataTP] = SimulateRxnKcatCoef(model_tmp,enzymedataSEC,enzymedataTP);  
    enzymedata_new = enzymedata;
    enzymedata_new.proteins = [enzymedata_new.proteins;enzymedataTP.proteins];
    enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs;enzymedataTP.proteinMWs];
    enzymedata_new.proteinLength = [enzymedata_new.proteinLength;enzymedataTP.proteinLength];
    enzymedata_new.proteinExtraMW = [enzymedata_new.proteinExtraMW;enzymedataTP.proteinExtraMW];
    enzymedata_new.kdeg = [enzymedata_new.kdeg;enzymedataTP.kdeg];%%
    enzymedata_new.proteinPST = [enzymedata_new.proteinPST;enzymedataTP.proteinPST];
    enzymedata_new.rxns = [enzymedata_new.rxns;enzymedataTP.rxns];
    enzymedata_new.rxnscoef = [enzymedata_new.rxnscoef;enzymedataTP.rxnscoef];
    enzymedata_new = CombineEnzymedata(enzymedata_new,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
    

    for j = 1:length(mulist)
        disp(num2str(j))
        mu = mulist(j);
        
        model_tmp = changeRxnBounds(model_tmp,'BIOMASS',mu,'b');
        name = [TP{i} '_',num2str(mu*100),'_PP'];
        fileName = writeLPGlc(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);% no extra constraint
        allname = [allname;{fileName}];
    end
end
writeclusterfileLP(allname,'sub_TP1')


%% reference condition to facilitate further identification of engenieering
% targets
% load model and param
load('enzymedata_PP.mat')
load('enzymedataMachine_PP.mat')
load('enzymedataSEC_PP.mat')
load('enzymedataDummyER_PP.mat')
load('pcSecPichia.mat')

mulist = [0.02:0.02:0.08,0.08:0.01:0.16,0.16:0.02:0.34];

model = setMediaPP(model,4);
% set carbon source
model = changeRxnBounds(model,'Ex_glc_D',-1000,'l');
model = changeRxnBounds(model,'BIOMASS',1000,'u');
model = changeRxnBounds(model,'LIPIDS',1000,'u');
model = changeRxnBounds(model,'PROTEINS',1000,'u');
model = changeRxnBounds(model,'STEROLS',1000,'u');

model = changeRxnBounds(model,'Ex_glyc',0,'l');
model = changeRxnBounds(model,'BIOMASS_glyc',0,'b');
model = changeRxnBounds(model,'LIPIDS_glyc',0,'b');
model = changeRxnBounds(model,'PROTEINS_glyc',0,'b');
model = changeRxnBounds(model,'STEROLS_glyc',0,'b');

model = changeRxnBounds(model,'Ex_meoh',0,'l');
model = changeRxnBounds(model,'BIOMASS_meoh',0,'b');
model = changeRxnBounds(model,'LIPIDS_meoh',0,'b');
model = changeRxnBounds(model,'PROTEINS_meoh',0,'b');
model = changeRxnBounds(model,'STEROLS_meoh',0,'b');

% set oxygen
model = changeRxnBounds(model,'Ex_o2',-1000,'l');

rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
model.ub(rxn) = 0;

% Set optimization
rxnID = 'Ex_glc_D'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.37; % g/gCDW, total protein content (0.37 for glucose condition, 0.40 for methanol condition)
% % tot_protein = 0.4; 

f_modeled_protein = extractModeledprotein(model,'BIOMASS','PROTEIN[c]'); %g/gProtein
% BIOMASS is pseudo_biomass_rxn_id in the GEM
% PROTEIN[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.04;

clear tot_protein f_modeled_protein;
factor_k = 1;

enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);

model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

for j = 1:length(mulist)
    disp(num2str(j))
    mu = mulist(j);
    
    model_tmp = changeRxnBounds(model,'BIOMASS',mu,'b');
    name = ['ref_',num2str(mu*100),'_PP'];
    fileName = writeLPGlc(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);% no extra constraint
    allname = [allname;{fileName}];
end
% write cluster file
writeclusterfileLP(allname,'sub_TPref')
% % cd ../


