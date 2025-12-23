%% simulation3CSource


%% PP

addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecPichia/')


% load model and param
load('enzymedata_PP.mat')
load('enzymedataMachine_PP.mat')
load('enzymedataSEC_PP.mat')
load('enzymedataDummyER_PP.mat')
load('pcSecPichia.mat')


mkdir('Simulate3CSource_Ppa');
cd Simulate3CSource_Ppa/;

%% Set model
% set medium
model = setMediaPP(model,2);

% set carbon source
% For each carbon source simulation, set the desired uptake reaction to a negative value 
% (e.g., -1000 for unlimited uptake) and set all other carbon sources to 0.
model = changeRxnBounds(model,'Ex_glc_D',-1000,'l');% glucose
model = changeRxnBounds(model,'BIOMASS',1000,'u');
model = changeRxnBounds(model,'Ex_glyc',0,'l');
model = changeRxnBounds(model,'BIOMASS_glyc',0,'b');
model = changeRxnBounds(model,'Ex_meoh',0,'l');
model = changeRxnBounds(model,'BIOMASS_meoh',0,'b');
% set oxygen
model = changeRxnBounds(model,'Ex_o2',-1000,'l');

%% Set optimization

% Modify 'rxnID' for each simulation depending on the chosen carbon source.
rxnID = 'Ex_glc_D'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.37; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'BIOMASS','PROTEIN[c]');
% BIOMASS is pseudo_biomass_rxn_id in the GEM
% PROTEIN[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;
factor_k = 1;
f_unmodelER = 0.040;
f_erm = 0.083;


%% Solve LPs
mu_list = [0.01:0.01:0.34];
allname = cell(0,1);
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

 allName = cell(0,1);
for i = 1:length(mu_list)
    mu = mu_list(i);
%     f_carbon = 5.244714732847007e-01*mu;
%     f_carbon = f_carbon/0.38067
    model_tmp = changeRxnBounds(model,'BIOMASS',mu,'b');
    disp(['mu = ' num2str(mu)]);
    name = [num2str(mu*100),'_GLC_PP'];
    fileName = writeLPGlc(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);
    allName{end+1} = fileName;
end

writeclusterfileLP(allName(:),'sub_1')

