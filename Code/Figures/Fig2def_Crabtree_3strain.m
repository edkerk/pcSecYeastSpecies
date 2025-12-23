%% Fig2def_Crabtree_3strain

%% KMX
addpath('../../Enzymedata/')
addpath('../../Model/')

% load model and param
load('enzymedata_KMX.mat')
load('enzymedataMachine_KMX.mat')
load('enzymedataSEC_KMX.mat')
load('enzymedataDummyER_KMX.mat')
load('pcSecKmarx.mat')

% Define growth rates and load experimental chemostat data
mu_list = [0.02:0.02:0.74];
fid = fopen('./ComplementaryData/chemostatDataKMX.tsv','r');
exp_data = textscan(fid,'%f32 %f32 %f32 %f32 %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4} exp_data{5}];
fclose(fid);
mod_data = zeros(numel(mu_list),5);
% Run simulations across dilution rates
for i = 1:length(mu_list)
    mu = num2str(mu_list(i)*100);
    fileName = strcat('Simulation_dilution',mu,'_GLC_KMX.lp.out');
    [sol_obj,sol_status,sol_full] = readSoplexResult(fileName,model);
    if isempty(sol_full) || all(sol_full == 0) || ~strcmp(sol_status, 'optimal')
        fprintf('Skipping invalid solution: %s (status: %s)\n', fileName, sol_status);
        continue;  
    end
    pos(1) = find(strcmp(model.rxnNames,'biomass exchange'));
    pos(2) = find(strcmp(model.rxnNames,'D-glucose exchange'));
    pos(3) = find(strcmp(model.rxnNames,'oxygen exchange'));
    pos(4) = find(strcmp(model.rxnNames,'carbon dioxide exchange'));
    pos(5) = find(strcmp(model.rxnNames,'ethanol exchange'));
    for j = 1:length(pos)
    mod_data(i,j) = abs(sol_full(pos(j)));
    end
end
% Save results
mod_data(mod_data == 0) = NaN;  
save('mod_data_KMX.mat',"mod_data");

%% SCE
% load model and param
load('enzymedata_SCE.mat')
load('enzymedataMachine_SCE.mat')
load('enzymedataSEC_SCE.mat')
load('enzymedataDummyER_SCE.mat')
load('pcSecYeast.mat')

% Define growth rates and load experimental chemostat data
mu_list = [0.02:0.02:0.44];
fid = fopen('./ComplementaryData/chemostatDataSCE.tsv','r');
exp_data = textscan(fid,'%f32 %f32 %f32 %f32 %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4} exp_data{5}];
fclose(fid);
mod_data = zeros(numel(mu_list),5);
% Run simulations across dilution rates
for i = 1:length(mu_list)
    mu = num2str(mu_list(i)*100);
    fileName = strcat('Simulation_dilution',mu,'_GLC_SCE.lp.out');
    [sol_obj,sol_status,sol_full] = readSoplexResult(fileName,model);
    if isempty(sol_full) || all(sol_full == 0) || ~strcmp(sol_status, 'optimal')
        fprintf('Skipping invalid solution: %s (status: %s)\n', fileName, sol_status);
        continue;  
    end
    pos(1) = find(strcmp(model.rxnNames,'growth'));
    pos(2) = find(strcmp(model.rxnNames,'D-glucose exchange'));
    pos(3) = find(strcmp(model.rxnNames,'oxygen exchange'));
    pos(4) = find(strcmp(model.rxnNames,'carbon dioxide exchange'));
    pos(5) = find(strcmp(model.rxnNames,'ethanol exchange'));
    for j = 1:length(pos)
    mod_data(i,j) = abs(sol_full(pos(j)));
    end
end
% Save results
mod_data(mod_data == 0) = NaN;  
save('mod_data_SCE.mat',"mod_data");


%% PP
% load model and param
load('enzymedata_PP.mat')
load('enzymedataMachine_PP.mat')
load('enzymedataSEC_PP.mat')
load('enzymedataDummyER_PP.mat')
load('pcSecPichia.mat')

% Define growth rates and load experimental chemostat data
mu_list = [0.01:0.01:0.33];
fid = fopen('./ComplementaryData/chemostatDataPP.tsv','r');
exp_data = textscan(fid,'%f32 %f32 %f32 %f32 %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4} exp_data{5}];
fclose(fid);
mod_data = zeros(numel(mu_list),5);

% Run simulations across dilution rates
for i = 1:length(mu_list)
    mu = num2str(mu_list(i)*100);
    fileName = strcat('Simulation_dilution',mu,'_GLC_PP.lp.out');
    [sol_obj,sol_status,sol_full] = readSoplexResult(fileName,model);
    if isempty(sol_full) || all(sol_full == 0) || ~strcmp(sol_status, 'optimal')
        fprintf('Skipping invalid solution: %s (status: %s)\n', fileName, sol_status);
        continue;  
    end
    pos(1) = find(strcmp(model.rxnNames,'Biomass exchange'));
    pos(2) = find(strcmp(model.rxnNames,'D-Glucose exchange'));
    pos(3) = find(strcmp(model.rxnNames,'Oxygen exchange'));
    pos(4) = find(strcmp(model.rxnNames,'CO2 exchange'));
    pos(5) = find(strcmp(model.rxnNames,'Ethanol exchange'));
    for j = 1:length(pos)
    mod_data(i,j) = abs(sol_full(pos(j)));
    end
end
% Save results
mod_data(mod_data == 0) = NaN;  
save('mod_data_PP.mat',"mod_data");


%% Figure
figure
hold on
set(gca, 'FontName', 'Arial', 'FontSize', 7,'FontWeight', 'bold');

% Plot model predictions (lines) and experimental data (markers)
cols = [0.16 0.44 0.56; 0.93 0.46 0.34; 0.60 0.80 0.20; 0.45 0.30 0.6];
b    = zeros(1,length(exp_data(1,:))-1);
for i = 1:length(exp_data(1,:))-1
    b(i) = plot(mod_data(:,1),mod_data(:,i+1),'Color',cols(i,:),'LineWidth',1.5);
    plot(exp_data(:,1),exp_data(:,i+1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerSize',3);
end
xlabel('Dilution rate [gDW/h]')
ylabel('Exchange fluxes [mmol/gDW/h]')

% Figure formatting
% % legend(b,'Glucose consumption','O2 consumption','CO2 production','ethanol production','Location','northwest','Box', 'off')
set(gcf,'units', 'centimeters','position',[10 10 5 5],'Color','none');
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1 1 3.5 3.5],'Color','none');
ymax = ylim;
max_val = ymax(2);
text(0.05, 0.9*max_val, 'Sce','FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
hold off

