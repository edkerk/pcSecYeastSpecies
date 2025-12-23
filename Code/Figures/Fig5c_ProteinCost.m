%% Fig5c_ProteinCost
% This script visualizes protein production costs in K. marxianus, 
% S. cerevisiae, and P. pastoris under glucose conditions.
addpath('../../Results/Protein_cost_TP/')

%% Load protein cost data
TP = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
tp_names = {'Gal';'GOX';'Inu';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'TRF';'Hb';'PHO';'AMS';'IP';'BGL'};
load('all_proteincost_gluKMX.mat')
glc_slope_finalKMX(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluSCE.mat')
glc_slope_finalSCE(:,:,1) = all_slope_glc(:,:,1);
load('all_proteincost_gluPP.mat')
glc_slope_finalPP(:,:,1) = all_slope_glc(:,:,1);


load("galactosidase_657_gluKMX.mat")
all_glc_slope_TPKMX(:,:,1) = res_slope_glc(a:b,:,1);
load("galactosidase_498_glu_SCE.mat")
all_glc_slope_TPSCE(:,:,1) = res_slope_glc(a:b,:,1);
load("PAS_chr1_galactosidase_352_glu_PP.mat")
all_glc_slope_TPPP(:,:,1) = res_slope_glc(a:b,:,1);

data = {abs(glc_slope_finalKMX(:,:,1)); abs(glc_slope_finalSCE(:,:,1)); abs(glc_slope_finalPP(:,:,1))};

tp_data = {abs(all_glc_slope_TPKMX(:,:,1));  
    abs(all_glc_slope_TPSCE(:,:,1));  
    abs(all_glc_slope_TPPP(:,:,1))    
};


%% Plot configuration
colors = [110,176,149; 200 139 168; 141,188,208]./255;  

tp_colors = [
    0.60 0.50 0.70; 0.45 0.55 0.80; 0.75 0.45 0.65; 
    0.50 0.70 0.60; 0.85 0.60 0.75; 0.40 0.60 0.75;
    0.70 0.65 0.50; 0.55 0.40 0.70; 0.80 0.70 0.40;
    0.65 0.75 0.55; 0.90 0.55 0.60; 0.50 0.50 0.60;
    0.70 0.50 0.55; 0.60 0.70 0.80; 0.30 0.40 0.60
];


figure;
hold on;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

for i = 1:length(data)
    current_data = data{i}(~isnan(data{i}));
    main_color = colors(i,:);
    box_color = min(main_color + 0.15, 1); 
    median_color = max(main_color - 0.15, 0); 
    outlier_color = min(main_color + 0.25, 1); 
     % Kernel density estimate for violin shape
    [f, yi] = ksdensity(current_data);
    f = f * 500;
    fill([f, -f(end:-1:1)]*0.5 + i, [yi, yi(end:-1:1)], ...
         main_color, 'FaceAlpha', 0.35, 'EdgeColor', main_color, 'LineWidth', 1);
    outli = isoutlier(current_data, 'quartiles');
    % Box plot overlay (quartiles, median, whiskers)
    qt25 = quantile(current_data, 0.25);
    qt75 = quantile(current_data, 0.75);
    med = median(current_data);
    plot([i i], [min(current_data(~outli)) qt25], 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    plot([i i], [qt75 max(current_data(~outli))], 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    fill(i + 0.1*[-1 1 1 -1], [qt25 qt25 qt75 qt75], ...
         box_color, 'FaceAlpha', 0.9, 'EdgeColor', max(main_color-0.1,0), 'LineWidth', 1);
    plot(i + 0.1*[-1 1], [med med], 'Color', median_color, 'LineWidth', 2.8);
end
%  Representative target proteins
all_tp_ys = [];
for j = 1:length(tp_names)
    for i = 1:length(tp_data)
        if j <= size(tp_data{i}, 1) && ~isnan(tp_data{i}(j,1))
            all_tp_ys(end+1) = tp_data{i}(j,1);
        end
    end
end
y_range = max(all_tp_ys) - min(all_tp_ys);
label_spacing = y_range * 0.12; 

tp_medians = zeros(length(tp_names),1);
for j = 1:length(tp_names)
    vals = [];
    for i = 1:length(tp_data)
        if j <= size(tp_data{i}, 1) && ~isnan(tp_data{i}(j,1))
            vals(end+1) = tp_data{i}(j,1);
        end
    end
    if ~isempty(vals)
        tp_medians(j) = median(vals);
    else
        tp_medians(j) = NaN;
    end
end
[~, sorted_idx] = sort(tp_medians);

used_ys = []; 
for j = sorted_idx'
    x_points = [];
    y_points = [];
    
    for i = 1:length(tp_data)
        if j <= size(tp_data{i}, 1) && ~isnan(tp_data{i}(j,1))
            x_points(end+1) = i;
            y_points(end+1) = tp_data{i}(j,1);
        end
    end
    
    if ~isempty(x_points)
        plot(x_points, y_points, '--', 'Color', [tp_colors(j,:) 0.7], ...
            'LineWidth', 1);
        
        scatter(x_points, y_points, 10, 'filled', 'MarkerEdgeColor', 'k', ...
               'MarkerFaceColor', tp_colors(j,:), 'LineWidth', 1);
        
        label_y = y_points(end)+200;
        if ~isempty(used_ys)
            while any(abs(label_y - used_ys) < label_spacing)
                label_y = label_y + label_spacing*0.7;
            end
        end
        used_ys(end+1) = label_y;
        
        end_x = 3.5;
        plot([x_points(end), end_x], [y_points(end), label_y], '--', ...
             'Color', [tp_colors(j,:) 0.4], 'LineWidth', 1);
        
        text(end_x, label_y, tp_names{j}, 'Color', tp_colors(j,:), ...
             'FontSize', 7, 'FontWeight', 'bold', ...
             'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    end
end
%% Axes and labels
set(gca,'units', 'centimeters', 'XTick', 1:length(data), 'XTickLabel', {'Kma'; 'Sce'; 'Ppa'},'LineWidth', 0.5, 'Position', [1.5 1 3.5 5]);
set(gcf, 'units', 'centimeters', 'position', [10 10 5.5 7]); 
ylabel({'Unit secretory cost','[mol glucose per mol protein]'}, 'FontSize', 7);
ylim([0 max(cellfun(@max, data)) * 1.2]);
xlim([0.5 4]);

