function model = changeGAM(model,GAM,NGAM)

bioPos = strcmp(model.rxnNames,'Biomass composition (g/g)');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','Phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

if nargin >2
    pos = strcmp(model.rxnNames,'ATP maintenance requirement');%NGAM
    model = setParam(model,'eq',model.rxns(pos),NGAM);% set both lb and ub
    model.lb(pos) = NGAM;
    model.ub(pos) = NGAM;
end

end