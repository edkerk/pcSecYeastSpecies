function [enzymedataSEC,modeled_ratio,meanprotein_info,missingsecP_ratio] = SimulateSecParam(model,protein_info,ProteinSequence)
% This function is to get all parameters for secretory pathway
% Input of this function is all protein abundance from literature (PMID 25346215)
% Use this to get the simulation for all parameters in the secretory
% pathway
% missingP_ratio is in percent of mg
% Yuesheng Zhang 20230908

% Load the proteome data
[num_abd, ~, raw_abd] = xlsread('protein_abundance_PP.xlsx','paxdb');
[num_mw, raw_mw, ~] = xlsread('protein_abundance_PP.xlsx','uniprot');

% Calculate total mass
abd_protein_list = raw_abd(2:end,2);
abd_abundance = raw_abd(2:end,3);
abd_abundance_mg = raw_abd(2:end,5);
% Load secretory complexes
[~,~,chap]=xlsread('TableS1_PP.xlsx','Secretory');
SecComplex = unique(chap(2:end,1),'stable');
Secgene = chap(2:end,1:2);
[~,Idx] = ismember(SecComplex,chap(:,1));
SecComplex_comp = chap(Idx,24);
SecComplex_coef = chap(Idx,25);
SecComplex_func = chap(Idx,26);

%% 1) get all rxnIDs for all proteins in proteome
rxnList = [];
geneError = [];
for i = 1:length(abd_protein_list)
    [model,allrxns,geneError_tmp] = addModificationforSec(model,abd_protein_list(i),protein_info,1);
    rxnList = [rxnList;allrxns];
    geneError = [geneError;geneError_tmp];
end

% match proteome abundance with rxnList
matchingList = cell(length(rxnList),3);

for i = 1:length(rxnList)
    [~, position] = find(rxnList{i,1} == '_', 3); 
    tmp = position(3);
    S{i,1} = extractBefore(rxnList{i}, tmp); 
    gene = {S{i,1}};

   [~,Idx] = ismember(gene,abd_protein_list);% This was changed [~,Idx] = ismember(gene,strrep(abd_protein_list,'-','_'))
    if Idx ~=0
        matchingList(i,:) = [gene,rxnList(i),abd_abundance(Idx)];
    else
        matchingList(i,:) = [gene,rxnList(i),0];
    end
end


% match the coef for the rxnList
enzymedata.proteins = setdiff(abd_protein_list,geneError);
enzymedata = calculateMW(enzymedata,ProteinSequence,protein_info);
enzymedata = getkdeg(enzymedata);
enzymedataSEC.enzyme = SecComplex;
enzymedataSEC.comp = SecComplex_comp;
enzymedataSEC.coefref = SecComplex_coef;
tmp = enzymedataSEC.coefref(strcmp(enzymedataSEC.enzyme,'sec_pdi1p_ero1p_complex')|strcmp(enzymedataSEC.enzyme,'sec_acc_Kar2p_complex')); % keep out accumulation rxn since it is not normal in cell
enzymedataSEC.coefref(strcmp(enzymedataSEC.enzyme,'sec_pdi1p_ero1p_complex')|strcmp(enzymedataSEC.enzyme,'sec_acc_Kar2p_complex')) = {'0*proteinLength'}; % keep out accumulation rxn since it is not normal in cell
a.rxns = rxnList; %  rxnlist for all proteins
enzymedata = SimulateRxnKcatCoef(a,enzymedataSEC,enzymedata);
[~,idx] = ismember(enzymedata.rxns,rxnList); % index the rxns with the sec complex
matchingList(idx(idx~=0),4) = num2cell(enzymedata.rxnscoef(idx~=0)); % coef
[~,idx] = ismember(enzymedata.rxns,rxnList); % index the rxns with the sec complex
matchingList(idx(idx~=0),4) = num2cell(enzymedata.rxnscoef(idx~=0)); % coef
clear enzymedata
enzymedataSEC.coefref(strcmp(enzymedataSEC.enzyme,'sec_pdi1p_ero1p_complex')|strcmp(enzymedataSEC.enzyme,'sec_acc_Kar2p_complex')) = tmp;


% taking care of the alternative pathway ERADM by HRD1 and DOA10, transloc
% by SEC61 or SSH
idx = find(contains(matchingList(:,2),'ERADM2'));
idx2 = find(contains(matchingList(:,2),'ERADM_sec_Ubc6p_Ubc7p_Hrd1p_Hrd3p_Der1p_complex'));%%changed in P.pastoris
matchingList(idx,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));
matchingList(idx2,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));
% SEC61/SSH1
idx = find(contains(matchingList(:,2),'_co_translation_TC_sec_SEC61C_complex'));
 idx2 = find(contains(matchingList(:,2),'_co_translation_TC_sec_SSH1C_complex'));
 matchingList(idx,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));
 matchingList(idx2,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));

% NG-golgi
idx = find(contains(matchingList(:,2),'_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pA_complex'));
idx2 = find(contains(matchingList(:,2),'_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pB_complex'));
idx3 = find(contains(matchingList(:,2),'_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pC_complex'));
matchingList(idx,3) = num2cell(cellfun(@(x) x/3,matchingList(idx,3)));
matchingList(idx2,3) = num2cell(cellfun(@(x) x/3,matchingList(idx,3)));
matchingList(idx3,3) = num2cell(cellfun(@(x) x/3,matchingList(idx,3)));


%coatother
 idx = find(contains(matchingList(:,2),'_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex'));
 idx2 = find(contains(matchingList(:,2),'_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Shl23p_Lst1p_Erv29p_Bet1p_Bos1p_complex'));
 matchingList(idx,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));
 matchingList(idx2,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));

%coatGPI
 idx = find(contains(matchingList(:,2),'_COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex'));
 idx2 = find(contains(matchingList(:,2),'_COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Shl23p_Lst1p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex'));
 matchingList(idx,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));
 matchingList(idx2,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));

%coatTransmembrane
 idx = find(contains(matchingList(:,2),'_COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Sec23p_Sec24p_Bet1p_Bos1p_complex'));
 idx2 = find(contains(matchingList(:,2),'_COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Shl23p_Lst1p_Bet1p_Bos1p_complex'));
 matchingList(idx,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));
 matchingList(idx2,3) = num2cell(cellfun(@(x) x/2,matchingList(idx,3)));

 %


idx_tmp = contains(model.rxns,'_complex_formation');
s_tmp = model.S(:,idx_tmp);
tf_tmp = s_tmp < 0;
max_subunit = max(sum(tf_tmp));

for i = 1:length(enzymedataSEC.enzyme)
    
    disp(['collect enzyme info' num2str(i) '/' num2str(length(enzymedataSEC.enzyme))]);
    
    enzyme_id = enzymedataSEC.enzyme{i};
    
    % add subunits
    enzfmtrxn_id = strcat(enzyme_id,'_formation');
    idx_tmp = ismember(model.rxns,enzfmtrxn_id);
    s_tmp = model.S(:,idx_tmp);
    subunits_tmp = model.mets(s_tmp < 0);
    na_tmp = repelem({''},max_subunit-length(subunits_tmp));
    subunits_tmp = cellfun(@(x) strrep(x,'_folding',''),subunits_tmp,'UniformOutput',false);
    subunits_tmp = cellfun(@(x) x(1:strfind(x,'[')-1),subunits_tmp,'UniformOutput',false);
    enzymedataSEC.subunit(i,:) = [subunits_tmp' na_tmp];
    
     % add stoichiometry of subunits
    stoichi_tmp = abs(full(s_tmp(s_tmp < 0)));
    na_tmp = repelem(0,max_subunit-length(stoichi_tmp));
    enzymedataSEC.subunit_stoichiometry(i,:) = [stoichi_tmp' na_tmp];
    
    % get subunit abun 
    [~,genes_index_proteome] = ismember(subunits_tmp,abd_protein_list);
    if any(genes_index_proteome ~= 0)
        for j = 1:length(subunits_tmp)
            if ~isequal(genes_index_proteome(j),0)    
                subunit_count(1,j) = abd_abundance(genes_index_proteome(j,1))';  
             elseif isequal(genes_index_proteome(j),0)
                 subunit_count(1,j) = num2cell(min(cell2mat(abd_abundance)));
            end
        end
    elseif genes_index_proteome == 0
        for j = 1:length(subunits_tmp)
        subunit_count(1,j) = num2cell(1000000/length(abd_protein_list));
        end
    end
    E0_tmp(i,1:length(subunits_tmp)) = subunit_count(1,1:length(subunits_tmp));   
    E0(i,1:length(subunits_tmp))= cell2mat(subunit_count(1,1:length(subunits_tmp)));
    enzymedataSEC.subunit_abun(i,1:length(subunits_tmp)) = subunit_count;
    % get catalyzed all enzyme abun * coef
    All_E = find(endsWith(matchingList(:,2),['_',enzyme_id]));
    E_sum(i,1:length(subunits_tmp)) = sum(cell2mat(matchingList(All_E,3)).*cell2mat(matchingList(All_E,4)));
    subunit_count = {};
end

%  ERAD should be only 30%  of total protein
u = 0.28; 
E_sum(strcmp(SecComplex_func,'ERAD'),:) = E_sum(strcmp(SecComplex_func,'ERAD'),:)*0.3;

% sum(V) <= Vsyn = kcat[E]
% (mu + kdeq)*sum([E]) <= kcat[E0]
allsecgene = setdiff(unique(enzymedataSEC.subunit(:)),'');
for i = 1:length(allsecgene)
    idx = ismember(enzymedataSEC.subunit,allsecgene(i));
    if length(find(idx)) > 1
        i
    E_sum(idx) = sum(E_sum(idx));
    end
end


kcat_tmp = (E_sum./E0).* enzymedataSEC.subunit_stoichiometry(:,1:length(E_sum(1,:)));

enzymedataSEC.kcat = median(kcat_tmp,2,'omitnan')*u*1.3; % 0.30 mean kdeg from the reference 
%enzymedataSEC.kcat = min(kcat_tmp,[],2)*(u+0.042);
%enzymedataSEC.proteins = strrep(setdiff(unique(enzymedataSEC.subunit(:)),''),'_','-'); % get all proteins involved in the sec pathway
enzymedataSEC.proteins = setdiff(unique(enzymedataSEC.subunit(:)),'');

% calculate modeled_protein coverage
list = endsWith(model.rxns,'_translation');
list = model.rxns(list);
list = strrep(list,'r_','');
modeled_proteins = strrep(list,'_peptide_translation','');

% calculate modeled protein ratio gram
[~,Idx] = ismember(modeled_proteins,abd_protein_list);
modeled_ratio = sum(cell2mat(abd_abundance_mg(Idx(Idx~=0))))/sum(cell2mat(abd_abundance_mg));%modeled_ratio = sum(abd_abundance_mg(Idx(Idx~=0)))/sum(abd_abundance_mg);

% proteins not take into account
missingprotein = setdiff(abd_protein_list,modeled_proteins);
[~,Idx] = ismember(missingprotein,abd_protein_list);
[~,Idx2] = ismember(missingprotein,protein_info(:,2));
missingprotein = missingprotein(Idx2~=0);
[~,Idx] = ismember(missingprotein,abd_protein_list);
[~,Idx2] = ismember(missingprotein,protein_info(:,2));

mean = sum(cell2mat(abd_abundance(Idx)).*cell2mat(protein_info(Idx2,3:9)))/sum(cell2mat(abd_abundance(Idx)));% cell2mat was added
mean_length = sum(cell2mat(abd_abundance(Idx)).*cell2mat(protein_info(Idx2,12)))/sum(cell2mat(abd_abundance(Idx)));
meanprotein_info = mean/mean_length*423;% 4.23 is the unmodeled protein secretory length*100 to save unnesaacy GTP and ATP

%% unmodeled protein secretory ratio
rxnList = [];
for i = 1:length(missingprotein)
    [model,allrxns] = addModificationforSec(model,missingprotein{i},protein_info,1);
    rxnList = [rxnList;allrxns];
end

% find unmodeled secretion protein list
for i = 1:length(rxnList)
    [~, position] = find(rxnList{i,1} == '_', 3); 
    tmp = position(3);
    S{i,1} = extractBefore(rxnList{i}, tmp);%added
    gene = {S{i,1}};
    [~,Idx] = ismember(gene,abd_protein_list);
    missingsecP(i) = gene;
end

% get missingSecP
missingsecP = unique(missingsecP);
[~,Idx] = ismember(missingsecP,abd_protein_list);
missingsecP_ratio = sum(cell2mat(abd_abundance_mg(Idx(Idx~=0))))/sum(cell2mat(abd_abundance_mg));


end
