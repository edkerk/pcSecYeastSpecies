function ecModel = tuningkcat(ecModel,enzymeName)
%tuning kcat
    prot_enzymeName = 'prot_'+string(enzymeName);
    rxns = findRxnsFromMets(ecModel, prot_enzymeName);
    ecModel.rxns_new = ecModel.rxns;
    ecModel.rxns = ecModel.ec.rxns;
    rxnIndices = findRxnIDs(ecModel, rxns);
    for i = 1:(length(rxnIndices)-1)
    ecModel.ec.kcat (rxnIndices(i,:)) = ecModel.ec.kcat (rxnIndices(i,:))*10;
    end
    ecModel.rxns_tunedkcat = ecModel.rxns;
    ecModel.rxns = ecModel.rxns_new;
end