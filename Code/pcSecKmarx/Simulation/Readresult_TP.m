%% Targets identification
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecKmarx/')


targetProteins = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};

res_list = cell('');
%% reference
% load model and param
load('enzymedata_KMX.mat')
load('enzymedataMachine_KMX.mat')
load('enzymedataSEC_KMX.mat')
load('enzymedataDummyER_KMX.mat')
load('pcSecKmarx.mat')

load('fluxesKMX_ref.mat')
mu = round(fluxes(ismember(model.rxns,'r_1913'),:),4);
fluxes = fluxes(:, (0.25 < mu) );
mu = round(fluxes(ismember(model.rxns,'r_1913'),:),4);
[a,b] = sort(mu,'ascend');
fluxes = fluxes(:,b);
nonzeroidx = any(fluxes);
fluxes = fluxes(:,nonzeroidx);
[abun_complex_id,abun_protein_id,abun_complex,abun_protein,~,~] = abundanceCalculation(model,fluxes);
mu = round(fluxes(ismember(model.rxns,'r_1913'),:),4);
idx_max = find(mu == max(mu));
idx_min = find(mu == min(mu));

res_ref.genelist = cell(0,1);
res_ref.r = zeros(0,1);
res_ref.p = zeros(0,1);
res_ref.slope = zeros(0,1);
res_ref.i = zeros(0,1);
res_ref.ratio = zeros(0,1);
res_ref.pro_abun = zeros(0,1);
for res_corted = 1:length(abun_protein_id)
    if ~any(abun_protein(res_corted,:) == 0)
        geneid = abun_protein_id(res_corted);
        res_ref.genelist = [res_ref.genelist;geneid];
        [RHOtmp,PVALtmp] = corr(abun_protein(res_corted,:)',mu','Type','Spearman');
        res_ref.r = [res_ref.r;RHOtmp];
        res_ref.p = [res_ref.p;PVALtmp];
        p = polyfit(abun_protein(res_corted,:),mu,1);
        res_ref.slope = [res_ref.slope;p(1)];
        res_ref.i = [res_ref.i;res_corted];
        x = abun_protein(res_corted,:);
        res_ref.pro_abun = [res_ref.pro_abun;x(idx_max)];
        res_ref.ratio = [res_ref.ratio;x(idx_max)/x(idx_min)];
    end
end
% remove ribosome
load('enzymedataMachine_KMX.mat')
ribosomegene = unique(enzymedataMachine.subunit(1:2,:));
ribosomegene = setdiff(ribosomegene,'');
[~,idx] = ismember(ribosomegene,res_ref.genelist);
res_ref.genelist(idx) = [];
res_ref.r(idx) = [];
res_ref.p(idx) = [];
res_ref.i(idx) = [];
res_ref.ratio(idx) = [];
res_ref.pro_abun(idx) = [];
res_ref.slope(idx) = [];
% label the genes
res_ref.priority(1:length(res_ref.r),1) = 0;
idx = find(res_ref.r > 0.9 & res_ref.ratio > 1.2);
res_ref.priority(idx) = 1; % similar with growth
idx = find(res_ref.r < -0.9 & res_ref.ratio < 0.86);
res_ref.priority(idx) = -1;% reverse with growth
% key note: recombiannt protein reverse with the growth, so identified
% targets for recombinant protein need to exclude the reverse targets for
% reference
%% target proteins
for i = 1:length(targetProteins)
    display([num2str(i) '/' num2str(length(targetProteins))]);
    
    load(['modelKMX_',targetProteins{i},'.mat']);
    load(['fluxesKMX_',targetProteins{i},'.mat']);
    tp_ex = [targetProteins{i},' exchange']; % index the tp production rate
    
    mu = round(fluxes(ismember(model.rxns,'r_1913'),:),4);
    fluxes = fluxes(:,mu > 0.25);
    mu = round(fluxes(ismember(model.rxns,'r_1913'),:),4);
    [a,b] = sort(mu,'ascend');
    fluxes = fluxes(:,b);
    mu = round(fluxes(ismember(model.rxns,'r_1913'),:),4);
    
    
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    q_tp = fluxes(ismember(model.rxns,tp_ex),:);
    fluxes = fluxes(:,find(q_tp));
    [abun_complex_id,abun_protein_id,abun_complex,abun_protein,~,~] = abundanceCalculation(model,fluxes);
    q_tp = fluxes(ismember(model.rxns,tp_ex),:);
    idx_max = find(q_tp == max(q_tp));
    idx_min = find(q_tp == min(q_tp));
    
    res.genelist = cell(0,1);
    res.r = zeros(0,1);
    res.p = zeros(0,1);
    res.slope = zeros(0,1);
    res.i = zeros(0,1);
    res.ratio = zeros(0,1);
    res.pro_abun = zeros(0,1);
    res.refstate = zeros(0,1);
    for res_corted = 1:length(abun_protein_id)
        if ~any(abun_protein(res_corted,:) == 0)
            geneid = abun_protein_id(res_corted);
            res.genelist = [res.genelist;geneid];
            [RHOtmp,PVALtmp] = corr(abun_protein(res_corted,:)',q_tp','Type','Spearman');
            res.r = [res.r;RHOtmp];
            res.p = [res.p;PVALtmp];
            p = polyfit(abun_protein(res_corted,:),q_tp,1);
            res.slope = [res.slope;p(1)];
            res.i = [res.i;res_corted];
            x = abun_protein(res_corted,:);
            res.pro_abun = [res.pro_abun;x(idx_max)];
            res.ratio = [res.ratio;x(idx_max)/x(idx_min)];
        end
    end
    
    % remove ribosome genes
    load('enzymedataMachine_KMX.mat')
    ribosomegene = unique(enzymedataMachine.subunit(1:2,:));
    ribosomegene = setdiff(ribosomegene,'');
    [~,idx] = ismember(ribosomegene,res.genelist);
    res.genelist(idx) = [];
    res.r(idx) = [];
    res.p(idx) = [];
    res.i(idx) = [];
    res.ratio(idx) = [];
    res.pro_abun(idx) = [];
    res.slope(idx) = [];
    % annotation and protein abundace ref from paxdb
    [num_abd, raw_abd, ~] = xlsread('protein_abundance_KM.xlsx','paxdb');
    raw_abd = raw_abd(2:end,:);
    res.ref_abun = zeros(length(res.genelist),1);
    [~,idx] = ismember(res.genelist,raw_abd(:,2));
    res.ref_abun(idx~=0) = num_abd(idx(idx~=0),6); % ref abundance in mmol/gDW
    
    [~, ~, raw_anno] = xlsread('Protein_annotation_uniprot.xlsx','uniprot');
    res.annot = cell(length(res.genelist),9);
    [~,idx] = ismember(res.genelist,raw_anno(:,1));
    res.annot(idx~=0,:) = raw_anno(idx(idx~=0),2:end);
    
    % annotation of homologous gene
    [~, ~, raw_homo] = xlsread('Yeast_homologue_genes_KM.xlsx');
    res.homo = cell(length(res.genelist),1);
    [~,idx] = ismember(res.genelist,raw_homo(:,2));
    res.homo(idx~=0,:) = raw_homo(idx(idx~=0),5);
    
    % annotation of protein wther it is metabolic
    res.met = repmat({'metabolic'},length(res.genelist),1);
    load('enzymedataSEC_KMX.mat')
    secprotein = unique(enzymedataSEC.subunit);
    [~,idx] = ismember(res.genelist,secprotein);
    res.met(idx ~= 0,1) = {'Secretory'};
        
    % difference between protein abundance/paxdb abundance
    res.difference = res.pro_abun./res.ref_abun;
    res.complex = nan;
    % check the state in reference condition without recombinant protein production
    [~,idx] = ismember(res.genelist,res_ref.genelist);
    res.refstate(idx~=0,1) = res_ref.priority(idx(idx ~=0));
    res.refstate(idx==0,1) = nan;
    
    res.priority = zeros(length(res.genelist),1);
    % priority rank for overexpression the higher, the better
    idx = find(res.r > 0.99 & res.ratio > 1);
    res.priority(idx) = 1;
    idx = find(res.r >= 0.99 & res.ratio > 1.2);
    res.priority(idx) = 2;
    idx = find(res.r >= 0.99 & res.ratio > 1.2 & res.difference >= 1);
    res.priority(idx) = 3;
    idx = find(res.r >= 0.99 & res.ratio > 1.2 & res.difference >= 1 & cellfun(@isempty,res.homo)& cellfun(@isempty,res.complex));
    res.priority(idx) = 4;
    idx = find(res.priority > 0 & res.refstate == -1); % those proteins may increse due to the reverse correlation with growth
    res.priority(idx) = 0.5;
    idx = find(res.priority > 0 & isnan(res.refstate)); % those proteins only occurs when recombinant protein production
    res.priority(idx) =  res.priority(idx) + 0.5;
    
    % priority rank for downregulation the lower, the better
    idx = find(res.r < -0.9 & res.ratio < 1);
    res.priority(idx) = -1;
    idx = find(res.r < -0.9 & res.ratio < 0.833);
    res.priority(idx) = -2;
    idx = find(res.r < -0.9 & res.ratio < 0.833 & res.difference < 0.5);
    res.priority(idx) = -3;
    idx = find(res.r < -0.9 & res.ratio < 0.833 & res.difference < 0.5 & cellfun(@isempty,res.homo)& cellfun(@isempty,res.complex));
    res.priority(idx) = -4;
    idx = find(res.priority < 0 & res.refstate == 1); % those proteins may decrease due to the correlation with growth
    res.priority(idx) = -0.5;
    idx = find(res.priority < 0 & isnan(res.refstate)); % those proteins are not present in the reference state, but with the change of recombinant protein
    res.priority(idx) =  res.priority(idx) - 0.5;
    
    
    t = table(res.genelist,res.priority,res.r,res.p,res.i,res.ratio,res.slope,res.pro_abun,res.ref_abun,res.difference,res.homo,res.met,res.complex,res.annot(:,1),res.annot(:,2),res.annot(:,3),res.annot(:,4),res.annot(:,5),res.annot(:,6),res.annot(:,7),res.annot(:,8),'VariableNames',{'protein' ,'priority','correlation', 'p_value','index_in_Model','ratio','slope','prot_abun(max_production)','ref_abun(paxdb)','abun_fold','homolougous_gene','metabolic/secretoy','complex/not',raw_anno{1,2:9}});
    res_list(1:length(res.genelist(res.priority >= 1)),i) = res.genelist(res.priority >= 1);
    res_listAnnotation(1:length(res.genelist(res.priority >= 1)),i) = res.annot(res.priority >= 1,3);
    res_listFunc(1:length(res.genelist(res.priority >= 1)),i) = res.met(res.priority >= 1);

    save(['res_',targetProteins{i},'_KMX_M2.mat'],'res');
    writetable(t,[targetProteins{i},'Target_KMX_M2.txt'],'Delimiter','\t','QuoteStrings',false,'WriteRowNames',true)
end
save('res_geneListKMX_M2.mat','res_list','targetProteins');


