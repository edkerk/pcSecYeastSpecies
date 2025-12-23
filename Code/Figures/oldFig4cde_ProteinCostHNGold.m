%% Fig4cde_ProteinCostHNG
% Compare unit secretory cost under NG-damaged conditions
TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'xynA';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
TP_abbr = {'Gal';'GOX';'Inu';'SOD';'Est';'HSA';'COL';'PH20';'xynA';'TRF';'Hb';'PHO';'AMS';'Insulin';'BGL'};
colors = [0.78 0.55 0.66; 0.65 0.75 0.83; 0.53 0.52 0.66];  
tp_colors = [0.6 0.6 0.6];

%% KMX
% Load protein cost data
load('all_proteincost_gluKMX.mat')
glc_slope_finalKMX(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluKMXHNG80.mat')
glc_slope_final_KMX80HNG(:,:,1) = all_slope_glc(:,:,1);
load('galactosidase_657_gluKMX.mat')
glc_slope_final_KMXTP(:,:,1) = res_slope_glc(a:b,:,1);
load('galactosidase_657_gluKMXHNG80.mat')
glc_slope_final_KMX80HNGTP(:,:,1) = res_slope_glc(a:b,:,1);

% Scatter plot with TPs 
x1 = abs(glc_slope_finalKMX(:,:,1));
y1 = abs(glc_slope_final_KMX80HNG(:,:,1));
threshold = 9.0;  
[x1, y1] = remove_outliers_mad(x1, y1, threshold);
x1_tp = abs(glc_slope_final_KMXTP(:,:,1));
y1_tp = abs(glc_slope_final_KMX80HNGTP(:,:,1));

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

h1 = scatter(x1, y1, 5, 'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', colors(1,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Kma');

h_tp = scatter(x1_tp, y1_tp, 10, 'MarkerFaceColor', tp_colors, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'TPs');

ylim([0,6000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory Cost(NG damaged)','[mol glucose per mol protein]'});


legend([h1, h_tp], {'Kma', 'TPs'}, 'Location', 'northwest', ...
    'Color', 'none', 'EdgeColor', 'none');


set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 
set(gca, 'units', 'centimeters', 'LineWidth', 1.5, 'Position', [1.5 1 3 4]);


%% SCE
% Load protein cost data
load('all_proteincost_gluSCE.mat')
glc_slope_finalSCE(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluSCEHNG80.mat')
glc_slope_final_SCE80HNG(:,:,1) = all_slope_glc(:,:,1);
load('galactosidase_498_gluSCE.mat')
glc_slope_final_SCETP(:,:,1) = res_slope_glc(a:b,:,1);
load('galactosidase_498_glu_SCEHNG80.mat')
glc_slope_final_SCE80HNGTP(:,:,1) = res_slope_glc(a:b,:,1);

% Scatter plot with TPs 
x2 = abs(glc_slope_finalSCE(:,:,1));
y2 = abs(glc_slope_final_SCE80HNG(:,:,1));
threshold = 9.0; 
[x2, y2] = remove_outliers_mad(x2, y2, threshold);
x2_tp = abs(glc_slope_final_SCETP(:,:,1));
y2_tp = abs(glc_slope_final_SCE80HNGTP(:,:,1));

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

h2 = scatter(x2, y2, 5, 'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', colors(2,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Sce');

h_tp = scatter(x2_tp, y2_tp, 10, 'MarkerFaceColor', tp_colors, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'TPs');


ylim([0,6000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory Cost(NG damaged)','[mol glucose per mol protein]'});

legend([h2, h_tp], {'Sce', 'TPs'}, 'Location', 'northwest', ...
    'Color', 'none', 'EdgeColor', 'none');
set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 
set(gca, 'units', 'centimeters', 'LineWidth', 1.5, 'Position', [1.5 1 3 4]);


%% PP
% Load protein cost data
load('all_proteincost_gluPP.mat')
glc_slope_finalPP(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluPPHNG80.mat')
glc_slope_final_PP80HNG(:,:,1) = all_slope_glc(:,:,1);
load('PAS_chr1_galactosidase_352_gluPP.mat')
glc_slope_final_PPTP(:,:,1) = res_slope_glc(a:b,:,1);
load('PAS_chr1_galactosidase_352_glu_PPHNG80.mat')
glc_slope_final_PP80HNGTP(:,:,1) = res_slope_glc(a:b,:,1);

% Scatter plot with TPs 
x3 = abs(glc_slope_finalPP(:,:,1));
y3 = abs(glc_slope_final_PP80HNG(:,:,1));
threshold = 9.0;  
[x3, y3] = remove_outliers_mad(x3, y3, threshold);
x3_tp = abs(glc_slope_final_PPTP(:,:,1));
y3_tp = abs(glc_slope_final_PP80HNGTP(:,:,1));

figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
h3 = scatter(x3, y3, 5, 'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', colors(3,:)*0.7, 'LineWidth', 1.5, 'DisplayName', 'Ppa');

h_tp = scatter(x3_tp, y3_tp, 10, 'MarkerFaceColor', tp_colors, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'TPs');


ylim([0,6000])
xlabel({'Unit secretory cost (native)','[mol glucose per mol protein]'});
ylabel({'Unit secretory Cost(NG damaged)','[mol glucose per mol protein]'});

legend([h3, h_tp], {'Ppa', 'TPs'}, 'Location', 'northwest', ...
    'Color', 'none', 'EdgeColor', 'none');
set(gcf, 'units', 'centimeters', 'position', [10 10 5 6]); 

set(gca, 'units', 'centimeters', 'LineWidth', 1.5, 'Position', [1.5 1 3 4]);


