function SimulateProteinCost(a,b)
% this function is to calculate the protein cost, ER cost of synthesizing each
% protein
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

% load model and param
load('enzymedata_KMX.mat')
load('enzymedataMachine_KMX.mat')
load('enzymedataSEC_KMX.mat')
load('enzymedataDummyER_KMX.mat')
load('pcSecKmarx.mat')

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
model = changeRxnBounds(model,'r_1878',0,'b');% pyruvate production

tot_protein = 0.45; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_1912','s_5006[c]'); %g/gProtein
% r_1912 is pseudo_biomass_rxn_id in the GEM
% s_5006[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;
factor_k = 1;
f_unmodelER = 0.046;
f_erm = 0.083;
allname = cell(0,1);
mulist = [0.05 0.1];
ratios = [5E-7,1E-6,5E-6,1E-5,2E-5];

for i = a:b
    i
    [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein(i),true);
    [enzymedataTP] = SimulateRxnKcatCoef(model_tmp,enzymedataSEC,enzymedataTP);
    enzymedata_new = enzymedata;
    enzymedata_new.proteins = [enzymedata_new.proteins;enzymedataTP.proteins];
    enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs;enzymedataTP.proteinMWs];
    enzymedata_new.proteinLength = [enzymedata_new.proteinLength;enzymedataTP.proteinLength];
    enzymedata_new.proteinExtraMW = [enzymedata_new.proteinExtraMW;enzymedataTP.proteinExtraMW];
    enzymedata_new.kdeg = [enzymedata_new.kdeg;enzymedataTP.kdeg];
    enzymedata_new.proteinPST = [enzymedata_new.proteinPST;enzymedataTP.proteinPST];
    enzymedata_new.rxns = [enzymedata_new.rxns;enzymedataTP.rxns];
    enzymedata_new.rxnscoef = [enzymedata_new.rxnscoef;enzymedataTP.rxnscoef];
    enzymedata_new = CombineEnzymedata(enzymedata_new,enzymedataSEC,enzymedataMachine,enzymedataDummyER);


    for k = 1:length(mulist)
        mu = mulist(k);
        model_tmp = changeRxnBounds(model_tmp,'r_1913',mu,'b');
        for j = 1:length(ratios)
            ratio = ratios(j);
            model_tmp = changeRxnBounds(model_tmp,[ERprotein{i},' exchange'],ratio,'b');
            rxnID = 'r_1726'; %minimize glucose uptake rate
            osenseStr = 'Maximize';
            name = [strrep(ERprotein{i},'new',''),'_',num2str(mu*100),'_ratio',num2str(1/ratio)];
            fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);
        
            allname = [allname;{fileName}];
        end
    end
end

% write cluster file
writeclusterfileLP(allname,['sub_proteinCost_',num2str(a)])

display([num2str(a) '/' num2str(length(ERprotein))]);
cd ../
end


