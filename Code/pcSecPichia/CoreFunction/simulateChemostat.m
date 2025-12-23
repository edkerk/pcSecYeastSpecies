function mod_data = simulateChemostat(model,exp_data)

%Relevant positions:
pos(1) = find(strcmp(model.rxns,'BIOMASS'));
pos(2) = find(strcmp(model.rxns,'Ex_glc_D'));
pos(3) = find(strcmp(model.rxns,'Ex_o2'));
pos(4) = find(strcmp(model.rxns,'Ex_co2'));

%Simulate chemostats:
mod_data = zeros(size(exp_data));
for i = 1:length(exp_data(:,1))
    %Fix biomass and minimize methanol:
    model = changeRxnBounds(model,model.rxns(pos(1)),exp_data(i,1));
    model = changeRxnBounds(model,model.rxns(pos(2)),-10,'l');
    model = changeObjective(model,model.rxns(pos(2)));
    sol   = optimizeCbModel(model,'max');

    %Store relevant variables:
    mod_data(i,:) = abs(sol.x(pos)');
end
end