%% simulation3CSource

%% SEC
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecYeast/');

% load model and param
load('enzymedata_SCE.mat');
load('enzymedataSEC_SCE.mat');
load('enzymedataDummyER_SCE.mat');
load('enzymedataMachine_SCE.mat')
load('pcSecYeast.mat');


mkdir('Simulate3CSource_Sce');
cd Simulate3CSource_Sce/;

%% Set model
% set medium
model = setMedia(model,2);

% set carbon source
% For each carbon source simulation, set the desired uptake reaction to a negative value 
% (e.g., -1000 for unlimited uptake) and set all other carbon sources to 0.
model = changeRxnBounds(model,'r_1714',-1000,'l');% GLC
model = changeRxnBounds(model,'r_1709',0,'l');% FRU
model = changeRxnBounds(model,'r_2058',0,'l');% SUC
model = changeRxnBounds(model,'r_1931',0,'l');% MAL
model = changeRxnBounds(model,'r_1710',0,'l');% GAL


% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_1810',0,'b');% glycine production
model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production

%% Set optimization

% Modify 'rxnID' for each simulation depending on the chosen carbon source.
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.046;
% % f_unmodelER =  0.046;

clear tot_protein f_modeled_protein;
factor_k = 1;

%% Solve LPs
mu_list = [0.01:0.01:0.44];
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

 allName = cell(0,1);
for i = 1:length(mu_list)
    mu = mu_list(i);
%     f_carbon = 5.244714732847007e-01*mu;
%     f_carbon = f_carbon/0.38067
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    name = [num2str(mu*100),'_GLC_SCE'];
    fileName = writeLPSCE(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);
    allName{i} = fileName;
end

writeclusterfileLP(allName(:),'sub_1')
