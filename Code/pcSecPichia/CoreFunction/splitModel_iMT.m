%% splitModel
function model_split = splitModel_iMT(model)
% split reactions with isozymes and reversible reactions used for general
% model

% Yu Chen (Nov 2021)
% Feiran Li (Nov 2021)

disp('Splitting reactions...');
if ~isfield(model,'grRules')
    model.grRules=rulesTogrrules(model);
end
index1 = find(ismember(model.grRules,'')); % no GPR

% a starting model only contains reactions that will not have enzyme constraints
model_split = struct();
model_split.rxns = model.rxns(index1);
model_split.mets = model.mets;
model_split.S = model.S(:,index1);
model_split.lb = model.lb(index1);
model_split.ub = model.ub(index1);
model_split.c = model.c(index1);
model_split.b = model.b;
model_split.comps = model.comps;
model_split.compNames = model.compNames;
model_split.rxnNames = model.rxnNames(index1);
%model_split.grRules = model.grRules(index1);%
model_split.genes = model.genes;
model_split.metNames = model.metNames;
model_split.metFormulas = model.metFormulas;
%model_split.metCharges = model.metCharges;
if ~isfield(model,'proteins')
    model_split.proteins=model.genes;
end
if isfield(model,'metMetaNetXID')
    model_split.metMetaNetXID = model.metMetaNetXID;
else
   model_split.metMetaNetXID = cell(length(model.mets),1); 
   model.metMetaNetXID = cell(length(model.mets),1);
end

if isfield(model,'modelID')
    model_split.modelID = [model.modelID,'(split)'];
else
    model_split.modelID = 'splited_model';
end

rxnidx = transpose(1:length(model.rxns));
rxnidx = setdiff(rxnidx,index1);

for m = 1:length(rxnidx)
    
    i = rxnidx(m);
    rxn = model.rxns{i};
    coeflist = full(model.S(:,i));
    substrates = model.mets(coeflist < 0);
    coef_sub = coeflist(coeflist < 0);
    products = model.mets(coeflist > 0);
    coef_pro = coeflist(coeflist > 0);
    lb = model.lb(i);
    ub = model.ub(i);
    reactionname_tmp = model.rxnNames{i};
    
    z = model.grRules{i};
    z = strsplit(z,' or ')';
    
    tfrvs = lb < 0;
    for j = 1:length(z)
        gr = z{j};
        gr = strrep(gr,'(','');gr = strrep(gr,')','');gr = strtrim(gr);
        if tfrvs % reversible
            % fwd
            if length(z) > 1
                rxnnew = [rxn,'_no_',num2str(j),'_fwd'];
            else
                rxnnew = [rxn,'_no_1_fwd'];
            end
            model_split = addReaction(model_split,rxnnew,...
                'reactionName',reactionname_tmp,...
                'metaboliteList',[substrates;products]',...
                'stoichCoeffList',[coef_sub;coef_pro]',...
                'lowerBound',0,'upperBound',ub,'geneRule',gr);
            % rvs
            if length(z) > 1
                rxnnew = [rxn,'_no_',num2str(j),'_rvs'];
            else
                rxnnew = [rxn,'_no_1_rvs'];
            end
            model_split = addReaction(model_split,rxnnew,...
                'reactionName',reactionname_tmp,...
                'metaboliteList',[substrates;products]',...
                'stoichCoeffList',[-coef_sub;-coef_pro]',...
                'lowerBound',0,'upperBound',-lb,'geneRule',gr);
        else
            if length(z) > 1
                rxnnew = [rxn,'_no_',num2str(j),'_fwd'];
            else
                rxnnew = [rxn,'_no_1_fwd'];
            end
            
            model_split = addReaction(model_split,rxnnew,...
                'reactionName',reactionname_tmp,...
                'metaboliteList',[substrates;products]',...
                'stoichCoeffList',[coef_sub;coef_pro]',...
                'lowerBound',0,'upperBound',ub,'geneRule',gr);
        end
    end
    
    if rem(m,1000) == 0 || m == length(rxnidx)
        disp(['-> Splitting reactions: ready with ' num2str(m) '/' num2str(length(rxnidx))]);
    end
end

% add proteins
[~,idx] = ismember(model_split.genes,model.genes);
% sort mets MNX IDs
[~,idx] = ismember(model_split.mets,model.mets);
model_split.metMetaNetXID = model.metMetaNetXID(idx);
[grRules,rxnGeneMat]   = standardizeGrRules(model_split,true);
model_split.rxnGeneMat = rxnGeneMat;
model_split.grRules = grRules;

end
