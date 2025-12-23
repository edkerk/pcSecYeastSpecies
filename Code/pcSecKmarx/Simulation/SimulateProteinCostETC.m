%% simulationProteinCostETC
% This function calculates protein cost and ER cost of protein synthesis under temperature-dependent conditions

% % KMX
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecKmarx/')

mkdir('SimulateProteinCost')
cd SimulateProteinCost

[~,~,protein_info] = xlsread('Protein_Information_KM.xlsx');
protein_info = protein_info(2:end,:);
ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strcat(ERprotein,'new');
protein_info(:,2) = strcat(protein_info(:,2),'new');
TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
ERprotein = [ERprotein;TP];

t_list = [30,40,50];

allName_ETC = cell(0,1);
for k = 1:length(t_list)
    t = t_list(k);
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
    tot_protein = 0.45;
    f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
    % r_1912 is pseudo_biomass_rxn_id in the GEM
    % s_5006[c] is protein id
    
    f = tot_protein * f_modeled_protein;
    f_unmodelER = 0.046;
    clear tot_protein f_modeled_protein;
    factor_k = 1;
    %% Solve LPs
    mulist = [0.05 0.1];
    ratios = [5E-7,1E-6,5E-6,1E-5,2E-5];
    model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

% %     set NGAM
    NGAM_T = getNGAMT(t);
    model = changeRxnBounds(model,'r_1905',NGAM_T,'b');

    for n = a:b
        n
        [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein(n),true);
        [enzymedataTP] = SimulateRxnKcatCoef(model_tmp,enzymedataSEC,enzymedataTP);
        enzymedata_new = enzymedataT;
        enzymedata_new.proteins = [enzymedata_new.proteins;enzymedataTP.proteins];
        enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs;enzymedataTP.proteinMWs];
        enzymedata_new.proteinLength = [enzymedata_new.proteinLength;enzymedataTP.proteinLength];
        enzymedata_new.proteinExtraMW = [enzymedata_new.proteinExtraMW;enzymedataTP.proteinExtraMW];
        enzymedata_new.kdeg = [enzymedata_new.kdeg;enzymedataTP.kdeg];
        enzymedata_new.proteinPST = [enzymedata_new.proteinPST;enzymedataTP.proteinPST];
        enzymedata_new.rxns = [enzymedata_new.rxns;enzymedataTP.rxns];
        enzymedata_new.rxnscoef = [enzymedata_new.rxnscoef;enzymedataTP.rxnscoef];
        enzymedata_new = CombineEnzymedata(enzymedata_new,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
        for m = 1:length(mulist)
            mu = mulist(m);
            model_tmp = changeRxnBounds(model_tmp,'r_1913',mu,'b');
            for j = 1:length(ratios)
                ratio = ratios(j);
                model_tmp = changeRxnBounds(model_tmp,[ERprotein{n},' exchange'],ratio,'b');            
                rxnID = 'r_1726'; %minimize glucose uptake rate   
                osenseStr = 'Maximize';
                name = [strrep(ERprotein{n},'new',''),'_',num2str(mu*100),'_ratio',num2str(1/ratio),'_',num2str(t)];
                fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);
                allName_ETC = [allName_ETC;{fileName}];
            end
        end
    end
end

writeclusterfileLP(allName_ETC(:),'sub_1')


