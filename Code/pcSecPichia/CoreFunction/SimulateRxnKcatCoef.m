function [enzymedata] = SimulateRxnKcatCoef(model,enzymedataSEC,enzymedata)
SecComplex = enzymedataSEC.enzyme;
SecComplex_coef = enzymedataSEC.coefref;
%% simulate coef for each rxn in protein modification
enzymedata.rxnscoef = [];
enzymedata.rxns = [];
for i = 1:length(SecComplex)
    rxns_by_this_complex = find(endsWith(model.rxns,['_',SecComplex{i}])& ~startsWith(model.rxns,'dummyER'));
    coefref_tmp = SecComplex_coef{i};
    %print kcat constriants in lp file
    rxnList = model.rxns(rxns_by_this_complex);
    %protList = extractBefore(rxnList,11);
    
    % fine gene with _A_ _B_
    protmp = {};% above was replaced by below
    for j = 1:length(rxnList)
        [~, position] = find(rxnList{j,1} == '_',3); 
        tmp = position(3);
        protmp{j,1} = extractBefore(rxnList{j}, tmp); 

    end
    protList = protmp;
    [~,Idx] = ismember(protList,enzymedata.proteins);
    if any(Idx)
        DSB = enzymedata.proteinPST(Idx(Idx~=0),strcmp(enzymedata.proteinPSTInfo,'DSB'));
        NG = enzymedata.proteinPST(Idx(Idx~=0),strcmp(enzymedata.proteinPSTInfo,'NG'));
        OG = enzymedata.proteinPST(Idx(Idx~=0),strcmp(enzymedata.proteinPSTInfo,'OG'));
        GPI = enzymedata.proteinPST(Idx(Idx~=0),strcmp(enzymedata.proteinPSTInfo,'GPI'));
        proteinLength = enzymedata.proteinLength(Idx(Idx~=0));
        if strcmp(coefref_tmp,'proteinLength')
            enzymedata.rxnscoef = [enzymedata.rxnscoef;enzymedata.proteinLength(Idx(Idx~=0))/467]; % bionumber 105224 average length
            enzymedata.rxns = [enzymedata.rxns;rxnList(Idx~=0)];
        elseif strcmp(coefref_tmp,'ProteinMW')
            enzymedata.rxnscoef = [enzymedata.rxnscoef;enzymedata.proteinMWs(Idx(Idx~=0))/54580]; % bionumber 115091 average weight
            enzymedata.rxns = [enzymedata.rxns;rxnList(Idx~=0)];
        else
            enzymedata.rxnscoef = [enzymedata.rxnscoef;eval(coefref_tmp)];
            enzymedata.rxns = [enzymedata.rxns;rxnList(Idx~=0)];
        end
    end
end
end