%% estimateModeledprotein 
function f_modeled_protein = estimateModeledprotein(model)

[num_abd, raw_abd, ~] = xlsread('protein_abundance_PP.xlsx','paxdb');
[num_mw, raw_mw, ~] = xlsread('protein_abundance_PP.xlsx','uniprot');

modeled_proteins = model.rxns(endsWith(model.rxns,'_translation'));
modeled_proteins = strrep(modeled_proteins,'r_','');
modeled_proteins = strrep(modeled_proteins,'_peptide_translation','');
%modeled_proteins = strrep(modeled_proteins,'_','-');

median_mw = median(num_mw); %used for proteins without molecular weight

% Calculate total mass
abd_protein_list = raw_abd(2:end,2);
abd_abundance = num_abd(:,3);

uniprot_protein_list = raw_mw(2:end,1);
uniprot_protein_mw = num_mw;

tot_mass = 0;
for i = 1:length(abd_protein_list)
    proteinid = abd_protein_list(i);
    if ismember(proteinid,uniprot_protein_list)
        mw = uniprot_protein_mw(ismember(uniprot_protein_list,proteinid));
    else
        mw = median_mw;
    end
    tot_mass = tot_mass + mw * abd_abundance(i);
end

modeled_mass = 0;
for i = 1:length(modeled_proteins)
    proteinid = modeled_proteins(i);
    if ismember(proteinid,uniprot_protein_list)
        mw = uniprot_protein_mw(ismember(uniprot_protein_list,proteinid));
    else
        i
        mw = median_mw;
    end
    
    if ismember(proteinid,abd_protein_list)
        abd = abd_abundance(ismember(abd_protein_list,proteinid));
    else
        abd = 0;
    end
    modeled_mass = modeled_mass + mw * abd;
end

f_modeled_protein = modeled_mass / tot_mass;
end
