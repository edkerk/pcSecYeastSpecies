%% changeBiomass 
function model = changeBiomass(model,f_modeled_protein,pseudo_biomass_rxn_id,protein_id)

% pseudo_biomass_rxn_id = 'r_1912';
% protein_id = 's_5006[c]';
rxnidx = ismember(model.rxns,pseudo_biomass_rxn_id);
% rxnidx = ismember(model.rxns,pseudo_protein_rxn_id);

% change protein coefficient


metidx = ismember(model.mets,protein_id);
% metidx = ismember(model.mets,carbohydrate_id);
% metidx = ismember(model.mets,lipid_id);
% metidx = ismember(model.mets,RNA_id);
% metidx = ismember(model.mets,DNA_id);

% metidx = ismember(model.mets,AA_id);

model.S(metidx,rxnidx) = -0.449998 * (1-f_modeled_protein);

