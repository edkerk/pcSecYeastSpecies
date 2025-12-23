%% simulationHNG

% % KMX
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecKmarx/')

mkdir('SimulateHNG')
cd SimulateHNG

allName_NG = cell(0,1);

load('pcSecKmarx_NG_KMX.mat','model');
load('enzymedata_NG_KMX.mat','enzymedata');
load('enzymedataSEC_NG_KMX.mat','enzymedataSEC');
load('enzymedataMachine_NG_KMX.mat','enzymedataMachine');
load('enzymedataDummyER_NG_KMX.mat','enzymedataDummyER');

[~,~,protein_info] = xlsread('protein_information_KM.xlsx');

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
model = setMediaKMX(model,2);
% set carbon source
model = changeRxnBounds(model,'r_1726',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1725',-1000,'l');
% block reactions
model = changeRxnBounds(model,'r_1759',0,'b');% acetate production
model = changeRxnBounds(model,'r_1758',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_1839',0,'b');% glycine production
model = changeRxnBounds(model,'r_1878',0,'b');% pyruvate production

%% Set optimization
rxnID = 'r_1726'; %minimize glucose uptake rate
osenseStr = 'Maximize';
tot_protein = 0.45;
f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
% r_1912 is pseudo_biomass_rxn_id in the GEM
% s_5006[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.046;
clear tot_protein f_modeled_protein;
factor_k = 1;

%% Solve LPs
mu_list = [0.01:0.02:0.7];
ratios = [0.1:0.1:1];
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

for m = 1:length(ratios)
    ratio = ratios(m);
    enzymedataNG_all = enzymedata_all;
    
    for i = 1:size(enzymedata_all.subunit, 1)
        subunit_genes = enzymedata_all.subunit(i, :);
        is_gene_present = any(ismember(subunit_genes, NG_protein));
    
        if is_gene_present
            enzymedataNG_all.kcat(i) = ratio * enzymedata_all.kcat(i);
        end
    end
      
    allName = cell(0,1);
    
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model_tmp = changeRxnBounds(model,'r_1913',mu,'b');
        disp(['mu = ' num2str(mu)]);
        name = ['NG_KMX',ratio,'_',num2str(mu*100)];
        fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedataNG_all,factor_k,name);    
        allName{end+1} = fileName;
    end
    allName_NG = [allName_NG, allName];
end

writeclusterfileLP(allName_NG(:),'sub_1')

