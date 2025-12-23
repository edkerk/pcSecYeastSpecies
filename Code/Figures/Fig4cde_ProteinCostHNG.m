%% Fig4cde_ProteinCostHNG
% Compare unit secretory cost under NG-damaged conditions
addpath('../../Results/Protein_cost_HNG/')

colors = [0.5294 0.7020 0.6510;   0.8902 0.5216 0.5804; 0.5490 0.7412 0.8196];  

%% Kma
% Load protein cost data
load('all_proteincost_gluKMX.mat')
glc_slope_finalKMX(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluKMXHNG80.mat')
glc_slope_final_KMX80HNG(:,:,1) = all_slope_glc(:,:,1);

x1 = abs(glc_slope_finalKMX(:,:,1));
y1 = abs(glc_slope_final_KMX80HNG(:,:,1));
threshold = 9.0;  
[x1, y1] = remove_outliers_mad(x1, y1, threshold);

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

h1 = scatter(x1, y1, 5, 'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', colors(1,:)*0.7, 'LineWidth', 1, 'DisplayName', 'Kma');

ylim([0,6000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory Cost(NG damaged)','[mol glucose per mol protein]'});


legend(h1, 'Kma', 'Location', 'northwest', ...
    'Color', 'none', 'EdgeColor', 'none');


set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 
set(gca, 'units', 'centimeters', 'LineWidth', 1, 'Position', [1.5 1 3 4]);


%% Sce
% Load protein cost data
load('all_proteincost_gluSCE.mat')
glc_slope_finalSCE(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluSCEHNG80.mat')
glc_slope_final_SCE80HNG(:,:,1) = all_slope_glc(:,:,1);

x2 = abs(glc_slope_finalSCE(:,:,1));
y2 = abs(glc_slope_final_SCE80HNG(:,:,1));
threshold = 9.0; 
[x2, y2] = remove_outliers_mad(x2, y2, threshold);

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

h2 = scatter(x2, y2, 5, 'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', colors(2,:)*0.7, 'LineWidth', 1, 'DisplayName', 'Sce');

ylim([0,6000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory Cost(NG damaged)','[mol glucose per mol protein]'});

legend(h2, 'Sce', 'Location', 'northwest', ...
    'Color', 'none', 'EdgeColor', 'none');
set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 
set(gca, 'units', 'centimeters', 'LineWidth', 1, 'Position', [1.5 1 3 4]);


%% Ppa
% Load protein cost data
load('all_proteincost_gluPP.mat')
glc_slope_finalPP(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluPPHNG80.mat')
glc_slope_final_PP80HNG(:,:,1) = all_slope_glc(:,:,1);
 
x3 = abs(glc_slope_finalPP(:,:,1));
y3 = abs(glc_slope_final_PP80HNG(:,:,1));
threshold = 9.0;  
[x3, y3] = remove_outliers_mad(x3, y3, threshold);


figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
h3 = scatter(x3, y3, 5, 'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', colors(3,:)*0.7, 'LineWidth', 1, 'DisplayName', 'Ppa');


ylim([0,6000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory Cost(NG damaged)','[mol glucose per mol protein]'});

legend(h3, 'Ppa', 'Location', 'northwest', ...
    'Color', 'none', 'EdgeColor', 'none');
set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 

set(gca, 'units', 'centimeters', 'LineWidth', 0.5, 'Position', [1.5 1 3 4]);


