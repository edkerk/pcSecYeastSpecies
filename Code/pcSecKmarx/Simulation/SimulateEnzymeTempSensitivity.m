%% SimulateEnzymeTempSensitivity

addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecKmarx/');

%% KMX

t_list = [20 25 30 35 40 45 50 55];
allName_ETC = cell(0,1);

for k = 1:length(t_list)
    t = t_list(k);

    load('enzymedata_KMX.mat')
    load('enzymedataMachine_KMX.mat')
    load('enzymedataSEC_KMX.mat')
    load('enzymedataDummyER_KMX.mat')
    load('pcSecKmarx.mat')
    load('ETC3mu.mat');
    idx = find(data.KM.t == t);        
    mu = data.KM.mu(idx);   
    
    [~,~,protein_info] = xlsread('protein_information_KM.xlsx');
    Tm_predKM = readtable('ProTstab2_resultTopt_KM.csv');
    Tm_predKM = table2struct(Tm_predKM, 'ToScalar', true);
    TmT90_infoKM = mapTmparams(Tm_predKM,enzymedata,protein_info);
    ETCparamsKM = get_ETCparams(TmT90_infoKM,t);

    gene_notmatch = {};
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
    
    [ETCparamsKM,enzymedataT] = calculate_kcatT(ETCparamsKM, enzymedata, 'KM',t);

    enzymedata_all = CombineEnzymedata(enzymedataT,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
    % %     set NGAM
    NGAM_T = getNGAMT(t);
    model = changeRxnBounds(model,'r_1905',NGAM_T,'b');
    
    allName = cell(0,1);
    % tuning kcatT
    for m = 1:length(enzymedata.proteins)
        ETCparamsKM_tuning = ETCparamsKM;
        enzymedataT_tuning = enzymedataT;
        genenumber = min(length(ETCparamsKM_tuning.Gene),length(enzymedataT_tuning.subunit));
        gene = ETCparamsKM_tuning.Gene(m);
        match_idx = find(strcmp(enzymedataT_tuning.subunit(:,1), gene));
        if ~isempty(match_idx)          
            for q = 1:length(match_idx)
                enzymedataT_tuning.kcat(match_idx(q),1) = 1.1*ETCparamsKM_tuning.kcatT(m);
            end
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
        rxnID = 'dilute_dummy'; %minimize glucose uptake rate
        osenseStr = 'Maximize';
        
        tot_protein = 0.45; %g/gCDW, estimated from the original GEM.
        f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
        % r_1912 is pseudo_biomass_rxn_id in the GEM
        % s_5006[c] is protein id
        
        f = tot_protein * f_modeled_protein;
        f_unmodelER = tot_protein * 0.046;
        % f_unmodelER = 0.0576;
        clear tot_protein f_modeled_protein;
        factor_k = 1;
        %% Solve LPs
        enzymedata_all_tuning = CombineEnzymedata(enzymedataT_tuning,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
        model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;
       
        model_tmp = changeRxnBounds(model,'r_1913',mu,'b');
        disp(['mu = ' num2str(mu)]);
        name = [num2str(t), '_',num2str(mu*100),'_',num2str(m),'_KMkcat1.1dummy'];
        fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all_tuning,factor_k,name);
        allName{end+1} = fileName;
    end
    model_tmp = changeRxnBounds(model,'r_1913',mu,'b');
    disp(['mu = ' num2str(mu)]);
    name = ['ref',num2str(t), '_',num2str(mu*100),'_KMkcatdummy'];
    fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);
    allName{end+1} = fileName;
    allName_ETC = [allName_ETC, allName];
end

writeclusterfileLP(allName_ETC(:),'sub_1')

