%% simulation3CSource

%% KMX
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecKmarx/')

% load model and param
load('enzymedata_KMX.mat')
load('enzymedataMachine_KMX.mat')
load('enzymedataSEC_KMX.mat')
load('enzymedataDummyER_KMX.mat')
load('pcSecKmarx.mat')

mkdir('Simulate3CSource_Kma');
cd Simulate3CSource_Kma/;

%% Set model
% set medium
model = setMediaKMX(model,2);

% set carbon source
% For each carbon source simulation, set the desired uptake reaction to a negative value 
% (e.g., -1000 for unlimited uptake) and set all other carbon sources to 0.
model = changeRxnBounds(model,'r_1726',-1000,'l');% GLC
model = changeRxnBounds(model,'r_1784',0,'l');% FRU
model = changeRxnBounds(model,'r_1890',0,'l');% SUC
model = changeRxnBounds(model,'r_1856',0,'l');% LAC
model = changeRxnBounds(model,'r_1785',0,'l');% GAL
model = changeRxnBounds(model,'r_1821',0,'l');% Inulin
model = changeRxnBounds(model,'r_1795',0,'l');% Xylose

% set oxygen
model = changeRxnBounds(model,'r_1725',-1000,'l');
% block reactions
model = changeRxnBounds(model,'r_1759',0,'b');% acetate production
model = changeRxnBounds(model,'r_1758',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_1878',0,'b');% pyruvate production


%% Set optimization

% Modify 'rxnID' for each simulation depending on the chosen carbon source.
rxnID = 'r_1726'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.45; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
% r_1912 is pseudo_biomass_rxn_id in the GEM
% s_5006[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = 0.046;
clear tot_protein f_modeled_protein;
factor_k = 1;

%% Solve LPs
mu_list = [0.01:0.01:0.75];
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

 allName = cell(0,1);
for i = 1:length(mu_list)
    mu = mu_list(i);
    model_tmp = changeRxnBounds(model,'r_1913',mu,'b');
    disp(['mu = ' num2str(mu)]);
    name = [num2str(mu*100),'_GLC_KMX'];
    fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);
    allName{i} = fileName;
end

writeclusterfileLP(allName(:),'sub_1')

