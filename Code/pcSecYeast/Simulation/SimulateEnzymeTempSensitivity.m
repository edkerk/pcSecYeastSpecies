%% SimulateEnzymeTempSensitivity

addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecYeast/');

%% SCE

t_list = [20 25 30 35 40 45];
allName_ETC = cell(0,1);

for k = 1:length(t_list)
    t = t_list(k);
    % load model
    load('enzymedata_SCE.mat')
    load('enzymedataMachine_SCE.mat')
    load('enzymedataSEC_SCE.mat')
    load('enzymedataDummyER_SCE.mat')
    load('pcSecYeast.mat')
    load('ETC3mu.mat');
    idx = find(data.SCE.t == t);        
    mu = data.SCE.mu(idx);   


    [~,~,protein_info] = xlsread('Protein_Information_SCE.xlsx');
    Tm_predSCE = readtable('ProTstab2_resultTopt_SCE.csv');
    Tm_predSCE = table2struct(Tm_predSCE, 'ToScalar', true);
    TmT90_infoSCE = mapTmparams(Tm_predSCE,enzymedata,protein_info);
    ETCparamsSCE = get_ETCparams(TmT90_infoSCE,t);
    
    gene_notmatch = {};
    for i = 1:length(enzymedata.proteins)
        gene_idx = find(strcmp(enzymedata.proteins{i}, ETCparamsSCE.Gene), 1);
        
        if isempty(gene_idx)
            gene_notmatch{end+1,1} = enzymedata.proteins{i};
             average_fUT = median(ETCparamsSCE.fUT(:));
            enzymedata.kdeg(i) = max(average_fUT, 0.3);
        else
            enzymedata.kdeg(i) = ETCparamsSCE.fUT(gene_idx);

        end
    end

    [ETCparamsSCE,enzymedataT] = calculate_kcatT(ETCparamsSCE, enzymedata, 'SCE',t);



    enzymedata_all = CombineEnzymedata(enzymedataT,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
% %      set NGAM
    NGAM_T = getNGAMT(t);
    model = changeRxnBounds(model,'r_4046',NGAM_T,'b');
   
    allName = cell(0,1);


    % tuning kcatT
    for m = 1:length(enzymedata.proteins)
        ETCparamsSCE_tuning = ETCparamsSCE;
        enzymedataT_tuning = enzymedataT;
    
        genenumber = min(length(ETCparamsSCE_tuning.Gene),length(enzymedataT_tuning.subunit));
        gene = ETCparamsSCE_tuning.Gene(m);
        match_idx = find(strcmp(enzymedataT_tuning.subunit(:,1), gene));
        if ~isempty(match_idx)          
            for q = 1:length(match_idx)
                enzymedataT_tuning.kcat(match_idx(q),1) = 1.1*ETCparamsSCE_tuning.kcatT(m);
            end
        end
    
        %% Set model
        % set medium
        model = setMedia(model,2);
        % set carbon source
        model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
        % set oxygen
        model = changeRxnBounds(model,'r_1992',-1000,'l');
        % block reactions
        model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
        model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production
        model = changeRxnBounds(model,'r_1810',0,'b');% glycine production
        model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production
        %% Set optimization
        rxnID = 'dilute_dummy'; %minimize glucose uptake rate
        osenseStr = 'Maximize';
        
        tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
        f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
        % r_4041 is pseudo_biomass_rxn_id in the GEM
        % s_3717[c] is protein id
        
        f = tot_protein * f_modeled_protein;
        f_unmodelER = tot_protein * 0.046;
        
        clear tot_protein f_modeled_protein;
        factor_k = 1;
        %% Solve LPs
        enzymedata_all_tuning = CombineEnzymedata(enzymedataT_tuning,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
        model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;
        model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
        disp(['mu = ' num2str(mu)]);
        name = [num2str(t), '_',num2str(mu*100),'_',num2str(m),'_SCEkcat1.1dummy'];
        fileName = writeLPref(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all_tuning,factor_k,name);
        allName{end+1} = fileName;\
    end
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    name = ['ref',num2str(t), '_',num2str(mu*100),'_SCEkcatdummy'];
    fileName = writeLPSCE(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);
    allName{end+1} = fileName;
    allName_ETC = [allName_ETC, allName];
end

writeclusterfileLP(allName_ETC(:),'sub_ref')

