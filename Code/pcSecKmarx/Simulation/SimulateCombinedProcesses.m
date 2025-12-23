%% SimulateCombinedProcesses

%% KMX
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecKmarx/')

mkdir('SimulateETC_KM')
cd('SimulateETC_KM')

t_list = [10 15 20 25 30 35 40 45 50 55 60];


%% Define experimental setups
use_kdegT = true;   % Temperature-dependent kdegT (default = 0.3)
use_kcatT = true;   % Temperature-dependent kcatT (default = enzymedata)
use_NGAMT = true;   % Temperature-dependent NGAM (default = fixed NGAM)

prefixes = {'KM_kcatT', 'KM_kdegT', 'KM_NGAM'};
settings = [use_kcatT, use_kdegT, use_NGAMT];

selectedPrefix = prefixes{settings};
if use_kdegT && use_kcatT && use_NGAMT
    prefix = 'KM_PCSEC';
end
allName_ETC = cell(0,1);
for k = 1:length(t_list)
    t = t_list(k);
    % load model
    load('enzymedata_KMX.mat')
    load('enzymedataMachine_KMX.mat')
    load('enzymedataSEC_KMX.mat')
    load('enzymedataDummyER_KMX.mat')
    load('pcSecKmarx.mat')

    [~,~,protein_info] = xlsread('protein_information_KM.xlsx');
    Tm_predKM = readtable('ProTstab2_resultTopt_KM.csv');
    Tm_predKM = table2struct(Tm_predKM, 'ToScalar', true);
    TmT90_infoKM = mapTmparams(Tm_predKM,enzymedata,protein_info);
    ETCparamsKM = get_ETCparams(TmT90_infoKM,t);
  
    gene_notmatch = {};
    if use_kdegT
        for i = 1:length(enzymedata.proteins)
            gene_idx = find(strcmp(enzymedata.proteins{i}, ETCparamsKM.Gene), 1);
            
            if isempty(gene_idx)
                gene_notmatch{end+1,1} = enzymedata.proteins{i};
                average_fUT = median(ETCparamsKM.fUT(:));
                enzymedata.kdeg(i) = max(average_fUT, 0.3);
            else
                enzymedata.kdeg(i) = ETCparamsKM.fUT(gene_idx);
            end
        end
    else
        enzymedata.kdeg(:) = 0.3;
    end
    
    if use_kcatT
        [ETCparamsKM,enzymedataT] = calculate_kcatT(ETCparamsKM, enzymedata, 'KM',t);
        enzymedata_all = CombineEnzymedata(enzymedataT,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
    else
        enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
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
% %     rxnID = 'dilute_dummy';
    osenseStr = 'Maximize';
   
    tot_protein = 0.45;%g/gCDW, estimated from the original GEM.
    f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
    % r_1912 is pseudo_biomass_rxn_id in the GEM
    % s_5006[c] is protein id
    
    f = tot_protein * f_modeled_protein;
    f_unmodelER = tot_protein * 0.046;
    clear tot_protein f_modeled_protein;
    factor_k = 1;
    %% Solve LPs
    mu_list = [0.02:0.02:0.74];

    model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

%    set NGAM
    if use_NGAMT
        NGAM_T = getNGAMT(t);
        model = changeRxnBounds(model,'r_1905',NGAM_T,'b');
    end

    allName = cell(0,1);
    
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model_tmp = changeRxnBounds(model,'r_1913',mu,'b');
        disp(['mu = ' num2str(mu)]);
        name = [prefix,'_',num2str(t), '_',num2str(mu*100)];
        fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);
        allName{i} = fileName;
    end
    allName_ETC = [allName_ETC, allName];

end
% 
writeclusterfileLP(allName_ETC(:),'sub_1')
