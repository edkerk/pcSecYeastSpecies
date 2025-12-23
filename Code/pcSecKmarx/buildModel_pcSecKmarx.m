%% buildModel
% Timing: ~  s

% Before run this script, please 1) check the annotation file is correct
% for all sce proteins
addpath('../../Data/pcSecKmarx/')


%load protein information
[~,~,protein_info] = xlsread('protein_information_KM.xlsx');
 
tic;
%% Import K.M GEM
% cd ../model/;
% org_model = readCbModel('Kluyveromyces_marxianus-GEM.xml');
% cd ../../;
% save('Kluyveromyces_marxianus-GEM.mat','org_model');
% load('K.M_GEM.mat');
load('K.M_GEM_changeformula1912.mat');

model = renameAllMets(model);

%% Modify the original model
[model_updated,energyResults,redoxResults] = modifyKMGEM(model);
%clear org_model;
model = model_updated;
%% Reformulate the orginal model
%   1. Reactions with isozymes (i.e., 'OR' case in GPR) will be copied,
%   each isozyme will be assigned to each copy. 
%   2. Reversible reactions will be split into forward and reverse ones.

model_split = splitModel(model);
model = model_split;

clear model_updated model_split;

%% reformulate the Metabolic part
load('Protein_Sequence.mat'); % this refers all protein id and sequence in S.ce
% Add complex formation reactions based on protein stoichiometry
% promiscuous = findPromiscuous(model_split);
load('Protein_stoichiometry.mat');% obtained from pdbe see Supplementary Methods for detail info
model = addComplexRxns(model,Protein_stoichiometry,protein_info,ProteinSequence);
model = renameAllMets(model);

%% reformulate the Sec part
% Add Sec complexes
[~,~,proteins]=xlsread('TableS1_KM.xlsx','Secretory');
model = addMachineryComplex(model,proteins,protein_info,Protein_stoichiometry,ProteinSequence);

%% reformulate the ribosome/assembly factor/proteasome part
[~,~,protein_machinery]=xlsread('TableS1_KM.xlsx','Machinery');
model = addMachineryComplex(model,protein_machinery,protein_info,Protein_stoichiometry,ProteinSequence);
model = renameAllMets(model);
% manually update complex formation reactions for some complexes.
% model = updateComplexformation(model,protein_info);

%% Collect kcats for enzymes
enzymedata = collectkcats(model);

% % % manually update kcats for some reactions.
% % enzymedata = updatekcats(enzymedata);
% kcat values were integrated directly using the GECKO3
% % load('enzymedata_KMX.mat')
list = endsWith(model.rxns,'_translation');
list = model.rxns(list);
list = strrep(list,'r_','');
enzymedata.proteins = strrep(list,'_peptide_translation',''); % all protein

% get kdeg info from petri's proteome data
enzymedata = getkdeg(enzymedata);
%enzymedata.kdeg(1:end) = 0.1;% using 0.1 for now

% Calculate molecular weight for each enzyme
enzymedata = calculateMW(enzymedata,ProteinSequence,protein_info);

% % % match kapp in the dataset % optinal
% % enzymedata = matchkappToKcat(enzymedata,kapp4);
% % 
% which is based on Sec protein abundance and the sec machinery abundance
[enzymedataSEC,modeled_ratio,meanprotein_info,missingsecP_ratio] = SimulateSecParam(model,protein_info,ProteinSequence);
enzymedataSEC = calculateMW(enzymedataSEC,ProteinSequence,protein_info);
enzymedata = SimulateRxnKcatCoef(model,enzymedataSEC,enzymedata);
% calculate the enzyme molecule weight of other machinery proteins such as
% ribosome and proteasome
enzymedataMachine = SimulateMachineParam(model);
enzymedataMachine = calculateMW(enzymedataMachine,ProteinSequence,protein_info);

%   % Change kcats extremely low for original enzymes no need to do this
%     lowkcat = 3600;
%     enzymedata.kcat(enzymedata.kcat< lowkcat) = lowkcat;
%% Add dummy complex reactions
% Dummy complex is assumed to be a part of metabolic protein pool.
% Note that the dummy complex is synthesized or diluted in the unit of
% mmol/gCDW/h.
[protein_info,~] = SimulateDummyERParam(model,meanprotein_info,protein_info,enzymedataSEC);
model = addDummyER(model,'r_1906','s_5006[c]',protein_info); 
[~,enzymedataDummyER] = SimulateDummyERParam(model,meanprotein_info,protein_info,enzymedataSEC);
% % enzymedataDummyER = calculateMW(enzymedataDummyER,ProteinSequence,protein_info);

% r_1906 is pseudo_protein_rxn_id in the GEM
% s_5006[c] is protein id
model = addDummy(model,'r_1906','s_5006[c]'); 
% close refolding and dilute for misfolding protein those reactions are for
% later usage
refold_list = model.rxns(contains(model.rxns,'refolding_'));
model = changeRxnBounds(model,refold_list,0,'b');
% misfold_dilute_list = model.rxns(contains(model.rxns,'dilution_misfolding'));
% model = changeRxnBounds(model,misfold_dilute_list,0,'b');

%% Change the original biomass equation
% Estimate modeled proteome
f_modeled_protein = estimateModeledprotein(model);

f_modeled_protein = f_modeled_protein + missingsecP_ratio; % small amount of protein has been added into ER
f_modeled_protein = floor(f_modeled_protein * 100) / 100;


% set medium
model = setMediaKMX(model,2);% 1 biomass 2 GAM

model = fitGAM(model);
% Change the biomass equation
model = changeBiomass(model,f_modeled_protein,'r_1912','s_5006[c]');
% r_1912 is pseudo_biomass_rxn_id in the original GEM
% s_5006[c] is protein id in the original GEM
% unmodeled_cofactor[c] is newly added metabolite for unmodeled cofactor

% change GAM and NGAM
GAM= 10;
NGAM = 1;
model = changeGAM(model,GAM,NGAM);

model = renameAllMets(model);

%% Save model
save('pcSecKmarx.mat','model');
save('enzymedata_KMX.mat','enzymedata');
save('enzymedataSEC_KMX.mat','enzymedataSEC');
save('enzymedataMachine_KMX.mat','enzymedataMachine');
save('enzymedataDummyER_KMX.mat','enzymedataDummyER');
%% Save model to Excel
model_excel = model;
model_excel.subSystems = cell(length(model_excel.rxns),1);
model_excel.mets = model_excel.metNames;
writeCbModel(model_excel,'xls','pcSecKmarx.xls');
clear model_excel;
toc;