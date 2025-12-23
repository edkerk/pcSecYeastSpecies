%% Fig3cde_ProteinCostTDM
% Compare protein secretory cost under temperature stress
addpath('../../Results/Protein_cost_TDM/')

%% Kma
% Load slopes (glucose cost per protein) at 40–52 °C
load('all_proteincost_gluKMX.mat')
glc_slope_final(:,:,1) = all_slope_glc(:,:,1);

load('all_proteincost_gluKMX40.mat')
glc_slope_final_40(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluKMX47.mat')
glc_slope_final_47(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluKMX52.mat')
glc_slope_final_52(:,:,1) = all_slope_glc(:,:,1);

x = abs(glc_slope_final(:,:,1));
y1 = abs(glc_slope_final_40(:,:,1));
y2 = abs(glc_slope_final_47(:,:,1));
y3 = abs(glc_slope_final_52(:,:,1));

colors = [0 0 0;  0.55 0.74 0.82;  0.89 0.52 0.58];  

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
threshold = 9.0; 
[x1, y1] = remove_outliers_mad(x, y1, threshold);
[x2, y2] = remove_outliers_mad(x, y2, threshold);
[x3, y3] = remove_outliers_mad(x, y3, threshold);

% Scatter plots and regression fits
h1 = scatter(x1, y1, 5, 'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', colors(1,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Topt (40°C)');
h2 = scatter(x2, y2, 5, 'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', colors(2,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Tmid (47°C)');
h3 = scatter(x3, y3, 5, 'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', colors(3,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Tmax (52°C)');


fit1 = fitlm(x1, y1);
fit2 = fitlm(x2, y2);
fit3 = fitlm(x3, y3);

plot(x1, fit1.Fitted, 'Color', colors(1,:)*0.9, 'LineWidth', 1,'LineStyle','--');
plot(x2, fit2.Fitted, 'Color', colors(2,:)*0.9, 'LineWidth', 1,'LineStyle','--');
plot(x3, fit3.Fitted, 'Color', colors(3,:)*0.9, 'LineWidth', 1,'LineStyle','--');

h_kmx = scatter(NaN, NaN, 1, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'DisplayName', 'Kma');

ylim([0,8000])

xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory cost','[mol glucose per mol protein]'});
legend([h_kmx,h1, h2, h3],'Location', 'northwest','Color', 'none','EdgeColor','none');
set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 

set(gca, 'units', 'centimeters', 'LineWidth', 1, 'Position', [1.5 1 3 4]);


%% Sce
% Load slopes at 30–40 °C
load('all_proteincost_gluSCE.mat')
glc_slope_final(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluSCE30.mat')
glc_slope_final_30(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluSCE35.mat')
glc_slope_final_35(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluSCE40.mat')
glc_slope_final_40(:,:,1) = all_slope_glc(:,:,1);

x = abs(glc_slope_final(:,:,1));
y1 = abs(glc_slope_final_30(:,:,1));
y2 = abs(glc_slope_final_35(:,:,1));
y3 = abs(glc_slope_final_40(:,:,1));

colors = [0 0 0;  0.55 0.74 0.82;  0.89 0.52 0.58];

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
threshold = 9.0;  
[x1, y1] = remove_outliers_mad(x, y1, threshold);
[x2, y2] = remove_outliers_mad(x, y2, threshold);
[x3, y3] = remove_outliers_mad(x, y3, threshold);

% Scatter plots and regression fits
h1 = scatter(x1, y1, 5, 'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', colors(1,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Topt (30°C)');
h2 = scatter(x2, y2, 5, 'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', colors(2,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Tmid (35°C)');
h3 = scatter(x3, y3, 5, 'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', colors(3,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Tmax (40°C)');

fit1 = fitlm(x1, y1);
fit2 = fitlm(x2, y2);
fit3 = fitlm(x3, y3);

plot(x1, fit1.Fitted, 'Color', colors(1,:)*0.9, 'LineWidth', 1,'LineStyle','--');
plot(x2, fit2.Fitted, 'Color', colors(2,:)*0.9, 'LineWidth', 1,'LineStyle','--');
plot(x3, fit3.Fitted, 'Color', colors(3,:)*0.9, 'LineWidth', 1,'LineStyle','--');

h_sce = scatter(NaN, NaN, 1, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'DisplayName', 'Sce');

ylim([0,8000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory cost','[mol glucose per mol protein]'});

legend([h_sce,h1, h2,h3],'Location', 'northwest','Color', 'none','EdgeColor','none');
set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]);

set(gca, 'units', 'centimeters', 'LineWidth', 1, 'Position', [1.5 1 3 4]);



%% Ppa
% Load slopes at 28–44 °C
load('all_proteincost_gluPP.mat')
glc_slope_final(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluPP28.mat')
glc_slope_final_30(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluPP36.mat')
glc_slope_final_40(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluPP44.mat')
glc_slope_final_44(:,:,1) = all_slope_glc(:,:,1);

x = abs(glc_slope_final(:,:,1));
y1 = abs(glc_slope_final_30(:,:,1));
y2 = abs(glc_slope_final_40(:,:,1));
y3 = abs(glc_slope_final_44(:,:,1));

colors = [0 0 0;  0.55 0.74 0.82;  0.89 0.52 0.58]; 

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
threshold = 9.0;  
[x1, y1] = remove_outliers_mad(x, y1, threshold);
[x2, y2] = remove_outliers_mad(x, y2, threshold);
[x3, y3] = remove_outliers_mad(x, y3, threshold);

% Scatter plots and regression fits
h1 = scatter(x1, y1, 5, 'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', colors(1,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Topt (28°C)');
h2 = scatter(x2, y2, 5, 'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', colors(2,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Tmid (36°C)');
h3 = scatter(x3, y3, 5, 'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', colors(3,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Tmax (44°C)');

fit1 = fitlm(x1, y1);
fit2 = fitlm(x2, y2);
fit3 = fitlm(x3, y3);

plot(x1, fit1.Fitted, 'Color', colors(1,:)*0.9, 'LineWidth', 1,'LineStyle','--');
plot(x2, fit2.Fitted, 'Color', colors(2,:)*0.9, 'LineWidth', 1,'LineStyle','--');
plot(x3, fit3.Fitted, 'Color', colors(3,:)*0.9, 'LineWidth', 1,'LineStyle','--');

h_pp = scatter(NaN, NaN, 1, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'DisplayName', 'Ppa');

ylim([0,8000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory cost','[mol glucose per mol protein]'});

legend([h_pp,h1, h2, h3],'Location', 'northwest','Color', 'none','EdgeColor','none');
set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 

set(gca, 'units', 'centimeters', 'LineWidth', 1, 'Position', [1.5 1 3 4]);


