%% SimulateTP

%% KMX
% max production rate under various growth rates on minimal media
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecKmarx/')


TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
mulist = [0.02:0.04:0.1,0.1:0.02:0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65];

mkdir('SimulateTP_KMX');
cd SimulateTP_KMX/;
allname = cell(0,1);

%% Solve LPs
for i = 1:length(TP)
    % load model and param
    load('enzymedata_KMX.mat')
    load('enzymedataMachine_KMX.mat')
    load('enzymedataSEC_KMX.mat')
    load('enzymedataDummyER_KMX.mat')
    load('pcSecKmarx.mat')

    %% Set model
    % set medium
    model = setMediaKMX(model,4);
    % set carbon source
    model = changeRxnBounds(model,'r_1726',-1000,'l');% glucose
    % set oxygen
    model = changeRxnBounds(model,'r_1725',-1000,'l');
    % block reactions
    model = changeRxnBounds(model,'r_1759',0,'b');% acetate production
    model = changeRxnBounds(model,'r_1758',0,'b');% acetaldehyde production
    model = changeRxnBounds(model,'r_1878',0,'b');% pyruvate production
    rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
    model.ub(rxn) = 0;
    %% Set optimization
    tot_protein = 0.45; %g/gCDW, estimated from the original GEM.
    f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
    % r_1912 is pseudo_biomass_rxn_id in the GEM
    % s_5006[c] is protein id
    
    f = tot_protein * f_modeled_protein;
    f_mito = 0.1;
    f_unmodelER = tot_protein * 0.046;
    f_erm = 0.083;
    clear tot_protein f_modeled_protein;
    
    factor_k = 1; % doesn't matter this is to tune the saturation
    rxnID = strcat(TP{i},' exchange');
% %     rxnID = 'dilute_dummy';
    osenseStr = 'Maximize';
    [model_tmp,enzymedataTP] = addTargetProtein(model,TP(i)); 
    % NOTE: For humanized glycosylation simulations, adjust protein MW in addTargetProtein()/calculateMW().
    model = model_tmp;
    save(['modelKMX_',TP{i},'.mat'],'model')
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
    
    for j = 1:length(mulist)
        disp(num2str(j))
        mu = mulist(j);
        model_tmp = changeRxnBounds(model_tmp,'r_1913',mu,'b');
        name = [TP{i} '_',num2str(mu*100),'_KMX'];
        fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);% no extra constraint
        allname = [allname;{fileName}];
    end
end
writeclusterfileLP(allname,'sub_TP1')

%% reference condition to facilitate further identification of engenieering
% targets
load('enzymedata_KMX.mat')
load('enzymedataMachine_KMX.mat')
load('enzymedataSEC_KMX.mat')
load('enzymedataDummyER_KMX.mat')
load('pcSecKmarx.mat')
GAM= 10;
NGAM = 1;
model = changeGAM(model,GAM,NGAM);

model = setMediaKMX(model,2);
% set carbon source
model = changeRxnBounds(model,'r_1726',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1725',-1000,'l');
% block reactions
model = changeRxnBounds(model,'r_1759',0,'b');% acetate production
model = changeRxnBounds(model,'r_1758',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_1878',0,'b');% pyruvate production

rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
model.ub(rxn) = 0;

% Set optimization
rxnID = 'r_1726'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.45; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
% r_1912 is pseudo_biomass_rxn_id in the GEM
% s_5006[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.046;

clear tot_protein f_modeled_protein;
factor_k = 1;

enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

for j = 1:length(mulist)
    disp(num2str(j))
    mu = mulist(j);
    model_tmp = changeRxnBounds(model,'r_1913',mu,'b');
    name = ['ref_',num2str(mu*100),'_KMX'];
    fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);% no extra constraint
    allname = [allname;{fileName}];
end
% write cluster file
writeclusterfileLP(allname,'sub_TPref')
% % cd ../


