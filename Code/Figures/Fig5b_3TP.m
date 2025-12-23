%% Fig5b_3TP
% Comparative protein production potential across three yeasts
% This script compares the maximum protein production rates of 15 recombinant targets
% in K. marxianus (Kma), S. cerevisiae (Sce), and K. phaffii (Ppa).
addpath('../../Results/Growth_rate_TP/')

% Define colors and target protein lists
color = [31,119,180; 214,39,40; 44,160,44; 148,103,189; 255,127,14; 227,119,194; 127,127,127; 188,189,34; 23,190,207; 140,86,75; 255,152,150; 174,199,232; 152,223,138; 196,156,148; 247,182,210]./255;
targetProteins = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
targetProteins_rb = {'Gal';'GOX';'Inu';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'TRF';'Hb';'PHO';'AMS';'IP';'BGL'};

q_tp_max_values        = zeros(length(targetProteins), 3); 
q_tp_max_values_ratio  = zeros(length(targetProteins), 3);
mu_tp                  = zeros(length(targetProteins), 3);

% Iterate over target proteins and extract flux distributions from models
for i = 1:length(targetProteins)
    display([num2str(i) '/' num2str(length(targetProteins))]);
    
    % Kma
    load(['modelKMX_', targetProteins{i}, '.mat']);
    load(['fluxesKMX_', targetProteins{i}, '.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:, nonzeroidx);
    tp_ex = [targetProteins{i}, ' exchange'];
    mu = round(fluxes(ismember(model.rxns, 'r_1913'), :), 2);
    [mu, b] = sort(mu, 'ascend');
    fluxes = fluxes(:, b);
    q_tp = fluxes(ismember(model.rxns, tp_ex), :);
    q_tp_max_values(i, 1) = max(q_tp); 
    mu_tp(i,1) = mu(q_tp == max(q_tp));
    q_tp_max_values_ratio(i,1) = q_tp_max_values(i, 1)/mu_tp(i,1);


    % Sce
    load(['modelSCE_', targetProteins{i}, '.mat']);
    load(['fluxesSCE_', targetProteins{i}, '.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:, nonzeroidx);
    tp_ex = [targetProteins{i}, ' exchange'];
    mu = round(fluxes(ismember(model.rxns, 'r_2111'), :), 2);
    [mu, b] = sort(mu, 'ascend');
    fluxes = fluxes(:, b);
    q_tp = fluxes(ismember(model.rxns, tp_ex), :);
    q_tp_max_values(i, 2) = max(q_tp); 
    mu_tp(i,2) = mu(q_tp == max(q_tp));
    q_tp_max_values_ratio(i,2) = q_tp_max_values(i, 2)/mu_tp(i,2);

    
    % Ppa
    load(['modelPP_', targetProteins{i}, '.mat']);
    load(['fluxesPP_', targetProteins{i}, '.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:, nonzeroidx);
    tp_ex = [targetProteins{i}, ' exchange'];
    mu = round(fluxes(ismember(model.rxns, 'BIOMASS'), :), 2);
    [mu, b] = sort(mu, 'ascend');
    fluxes = fluxes(:, b);
    q_tp = fluxes(ismember(model.rxns, tp_ex), :);
    q_tp_max_values(i, 3) = max(q_tp); 
    mu_tp(i,3) = mu(q_tp == max(q_tp));
    q_tp_max_values_ratio(i,3) = q_tp_max_values(i, 3)/mu_tp(i,3);


end

%% Visualization
figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

% Assign distinct colors for species
color_KMX = [110,176,149] / 255;  
color_SCE = [227,132,147] / 255;  
color_PP  = [141,188,208] / 255;  
% Grouped bar plot of normalized protein production rates
bar_data = bar(q_tp_max_values_ratio, 'grouped');
bar_data(1).FaceColor = color_KMX;  
bar_data(2).FaceColor = color_SCE; 
bar_data(3).FaceColor = color_PP;  
for k = 1:length(bar_data)
    bar_data(k).EdgeColor = 'none';   
    bar_data(k).BarWidth = 1;
end
% Axis and labels
xticks(1:length(targetProteins_rb));
xticklabels(targetProteins_rb);
ylabel('Protein production rate [mmol/gDW]');
ylim([0, max(q_tp_max_values_ratio(:))*1.3]);
% Legend and figure formatting
legend('Kma', 'Sce', 'Ppa', 'FontSize', 7, 'Box', 'off');
set(gcf, 'units', 'centimeters', 'position', [10 10 8 7]); 
set(gca, 'units', 'centimeters', 'LineWidth', 0.5, 'Position', [1.5 1.5 6 5]);


