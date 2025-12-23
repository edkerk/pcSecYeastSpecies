%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = fitGAM(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = fitGAMforGlc(model)
%Set NGAM according to the Literature
model = changeRxnBounds(model,'ATPM',2.81,'l');
model = changeRxnBounds(model,'Ex_glyc',0);
model = changeRxnBounds(model,'Ex_meoh',0);
model = changeRxnBounds(model,'Ex_glc_D',-10,'l');
model = changeRxnBounds(model,'LIPIDS_glyc',0);
model = changeRxnBounds(model,'PROTEINS_glyc',0);
model = changeRxnBounds(model,'STEROLS_glyc',0);
model = changeRxnBounds(model,'BIOMASS_glyc',0);
model = changeRxnBounds(model,'LIPIDS_meoh',0);
model = changeRxnBounds(model,'PROTEINS_meoh',0);
model = changeRxnBounds(model,'STEROLS_meoh',0);
model = changeRxnBounds(model,'BIOMASS',1000,'u');
model = changeRxnBounds(model,'LIPIDS',1000,'u');
model = changeRxnBounds(model,'PROTEINS',1000,'u');
model = changeRxnBounds(model,'STEROLS',1000,'u');

%Load chemostat data:
fid = fopen('chemostatDataPP.tsv','r');
exp_data_all = textscan(fid,repmat('%f',1,5),'Delimiter','\t','HeaderLines',1);
exp_data = [exp_data_all{1:4}];
fclose(fid);

%GAMs to span:
disp('Estimating GAM:')
GAM = 5:5:70;

%1st iteration:
disp('1st iteration')
GAM = iteration(model,GAM,exp_data);

%2nd iteration:
disp('2nd iteration')
GAM = iteration(model,GAM-10:1:GAM+10,exp_data);

%3rd iteration:
disp('3nd iteration')
GAM = iteration(model,GAM-1:0.1:GAM+1,exp_data);

model = changeGAM(model,GAM,2.81);

%Plot fit:
mod_data = simulateChemostat(model,exp_data);
figure
hold on
cols = [0,1,0;0,0,1;1,0,0];
b    = zeros(1,length(exp_data(1,:))-1);
for i = 1:length(exp_data(1,:))-1
    b(i) = plot(mod_data(:,1),mod_data(:,i+1),'Color',cols(i,:),'LineWidth',2);
    plot(exp_data(:,1),exp_data(:,i+1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:))
end
title('GAM fitting for growth on glucose')
xlabel('Dilution rate [1/h]')
ylabel('Exchange fluxes [mmol/gDWh]')
legend(b,'Glucose consumption','O2 consumption','CO2 production','Location','northwest')
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GAM = iteration(model,GAM,exp_data)

fitting = ones(size(GAM))*1000;

for i = 1:length(GAM)
    %Modify GAM:
    model_i = changeGAM(model,GAM(i),2.81);%NGAM=2.81 (Glucose)
    disp(['changed GAM:' num2str(GAM(i))])
    %Simulate model and calculate fitting:
    mod_data   = simulateChemostat(model_i,exp_data);
    R          = (mod_data - exp_data)./exp_data;
    fitting(i) = sqrt(sum(sum(R.^2)));
    disp(['GAM = ' num2str(GAM(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[~,best] = min(fitting);

if best == 1 || best == length(GAM)
    error('GAM found is sub-optimal: please expand GAM search bounds.')
else
    GAM = GAM(best);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = changeGAM(model,GAM,NGAM)

bioPos = strcmp(model.rxns,'BIOMASS');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'atp[c]','adp[c]','h2o[c]','h[c]','pi[c]'},model.mets{i})) == 1;
    if S_ix ~= 0 && isGAM == 1
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

if nargin >1
    model = changeRxnBounds(model,'ATPM',NGAM,'l');    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    model = changeRxnBounds(model,model.rxns(pos(1)),exp_data(i,1),'l');
    model = changeRxnBounds(model,model.rxns(pos(2)),-10,'l');
    model = changeObjective(model,model.rxns(pos(2)));
    sol   = optimizeCbModel(model,'max');

    %Store relevant variables:
    mod_data(i,:) = abs(sol.x(pos)');
end
end



