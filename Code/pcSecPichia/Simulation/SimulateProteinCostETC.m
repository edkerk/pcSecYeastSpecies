%% simulationProteinCostETC
% This function calculates protein cost and ER cost of protein synthesis under temperature-dependent conditions

% % PP
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecPichia/');

mkdir('SimulateProteinCost')
cd SimulateProteinCost

[~,~,protein_info] = xlsread('Protein_Information_PP.xlsx');
protein_info = protein_info(2:end,:);
ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strcat(ERprotein,'new');
protein_info(:,2) = strcat(protein_info(:,2),'new');
TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
TP = "PAS_chr1_" + TP;
ERprotein = [ERprotein;TP];
t_list = [30 40 44];

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
    
    [ETCparamsPP,enzymedataT] = calculate_kcatT(ETCparamsPP, enzymedata, 'PP',t);

            
    %% Set model
    % set medium
    model = setMediaPP(model,2);
    % set carbon source
    model = changeRxnBounds(model,'Ex_glc_D',-1000,'l');% glucose
    % set oxygen
    model = changeRxnBounds(model,'Ex_o2',-1000,'l');
    
    %% Set optimization
    tot_protein = 0.37;
    f_modeled_protein = extractModeledprotein(model,'BIOMASS','PROTEIN[c]'); %g/gProtein
    % BIOMASS is pseudo_biomass_rxn_id in the GEM
    % PROTEIN[c] is protein id
    
    f = tot_protein * f_modeled_protein;
    f_unmodelER = 0.04;
    clear tot_protein f_modeled_protein;
    factor_k = 1;
    %% Solve LPs
    mulist = [0.05 0.1];
    ratios = [5E-7,1E-6,5E-6,1E-5,2E-5];
    model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

% %     set NGAM
    NGAM_T = getNGAMT(t);
    model = changeRxnBounds(model,'ATPM',NGAM_T,'b');

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
            model_tmp = changeRxnBounds(model_tmp,'BIOMASS',mu,'b');
            for j = 1:length(ratios)
                ratio = ratios(j);
                model_tmp = changeRxnBounds(model_tmp,[ERprotein{n},' exchange'],ratio,'b');            
            % Set optimization
                rxnID = 'Ex_glc_D'; %minimize glucose uptake rate    
                osenseStr = 'Maximize';           
                name = [strrep(ERprotein{n},'new',''),'_',num2str(mu*100),'_ratio',num2str(1/ratio),'_',num2str(t)];
                fileName = writeLPGlc(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);
            
                allName_ETC = [allName_ETC;{fileName}];
            end
        end
    end
end

writeclusterfileLP(allName_ETC(:),'sub_1')

