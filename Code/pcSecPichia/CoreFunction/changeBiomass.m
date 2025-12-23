%% changeBiomass 
function model = changeBiomass(model,f_modeled_protein,pseudo_biomass_rxn_id,protein_id)


rxnidx = ismember(model.rxns,pseudo_biomass_rxn_id);

% change protein coefficient
metidx = ismember(model.mets,protein_id);
proteincoef = model.S(metidx,rxnidx);
model.S(metidx,rxnidx) = 1 * proteincoef * (1-f_modeled_protein);
end


