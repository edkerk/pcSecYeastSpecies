%% SimulateCombinedProcesses

%%  PP
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecPichia/');

t_list = [10 15 20 25 30 35 40 45];
%% Define experimental setups
use_kdegT = true;   % Temperature-dependent kdegT (default = 0.3)
use_kcatT = true;   % Temperature-dependent kcatT (default = enzymedata)
use_NGAMT = true;   % Temperature-dependent NGAM (default = fixed NGAM)

prefixes = {'PP_kcatT', 'PP_kdegT', 'PP_NGAM'};
settings = [use_kcatT, use_kdegT, use_NGAMT];

selectedPrefix = prefixes{settings};
if use_kdegT && use_kcatT && use_NGAMT
    prefix = 'PP_PCSEC';
end
allName_ETC = cell(0,1);
for k = 1:length(t_list)
    t = t_list(k);
    % load model
    load('enzymedata_PP.mat')
    load('enzymedataMachine_PP.mat')
    load('enzymedataSEC_PP.mat')
    load('enzymedataDummyER_PP.mat')
    load('pcSecPichia.mat')

    [~,~,protein_info] = xlsread('Protein_Information_PP.xlsx');
    Tm_predPP = readtable('ProTstab2_resultTopt_PP.csv');
    Tm_predPP = table2struct(Tm_predPP, 'ToScalar', true);
    
    TmT90_infoPP = mapTmparams(Tm_predPP,enzymedata,protein_info);
    ETCparamsPP = get_ETCparams(TmT90_infoPP,t);
  
    
    gene_notmatch = {};
    if use_kdegT
        for i = 1:length(enzymedata.proteins)
            gene_idx = find(strcmp(enzymedata.proteins{i}, ETCparamsPP.Gene), 1);
            
            if isempty(gene_idx)
                gene_notmatch{end+1,1} = enzymedata.proteins{i};
                average_fUT = median(ETCparamsPP.fUT(:));
                enzymedata.kdeg(i) = max(average_fUT, 0.3);
            else
                enzymedata.kdeg(i) = ETCparamsSCE.fUT(gene_idx);
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
    model = setMediaPP(model,2);
    % set carbon source
    model = changeRxnBounds(model,'Ex_glc_D',-1000,'l');% glucose
    % set oxygen
    model = changeRxnBounds(model,'Ex_o2',-1000,'l');
    rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
    model.ub(rxn) = 0;    
    %% Set optimization
    rxnID = 'Ex_glc_D'; %minimize glucose uptake rate
% %     rxnID = 'dilute_dummy';
    osenseStr = 'Maximize';
    
    tot_protein = 0.37; % g/gCDW, total protein content (0.37 for glucose condition, 0.40 for methanol condition)
% %     tot_protein = 0.4; 
    f_modeled_protein = extractModeledprotein(model,'BIOMASS','PROTEIN[c]'); %g/gProtein
    % BIOMASS is pseudo_biomass_rxn_id in the GEM
    % PROTEIN[c] is protein id
    
    f = tot_protein * f_modeled_protein;
    f_unmodelER = tot_protein * 0.04;
    clear tot_protein f_modeled_protein;
    factor_k = 1;
    %% Solve LPs
    mu_list = [0.02:0.02:0.34];
    model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

%   set NGAM
    if use_NGAMT
        NGAM_T = getNGAMT(t);
        model = changeRxnBounds(model,'r_1905',NGAM_T,'b');
    end
    allName = cell(0,1);
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model_tmp = changeRxnBounds(model,'BIOMASS',mu,'b');
        disp(['mu = ' num2str(mu)]);
        name = [prefix,'_',num2str(mu*100),'_',num2str(t)];
        fileName = writeLPGlc(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);
        allName{i} = fileName;
    end
    allName_ETC = [allName_ETC, allName];

end
% 
writeclusterfileLP(allName_ETC(:),'sub_PP')
