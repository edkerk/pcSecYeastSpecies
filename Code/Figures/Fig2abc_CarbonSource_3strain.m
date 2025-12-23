%% Fig2abc_CarbonSource_3strain.m
%  Growth rate comparison between experimental and predicted data
%  Yeast species: K. marxianus (Kma), S. cerevisiae (Sce), P. pastoris (Ppa)
addpath('../../Results/CSource/')

legend_handles = [];
legend_names = [];
%% Kma
% Load table and remove rows with missing values in experimental/predicted rates
T = readtable('CSource_res.xlsx', 'Sheet', 'Kma');
rows_with_content = ~ismissing(T{:,6});
T(rows_with_content, :) = [];
carbon_sources = T.("Carbon_source");
mu_exp = T.("mu_in_vivo_h_1_");
mu_pred = T.("mu_in_silico_h_1_");

% Define plotting styles (color, marker, label) for carbon sources
styles = struct( ...
    'glucose',   struct('Color',[0, 0.45, 0.74],'Marker','o','Name','GLC'), ...
    'fructose',  struct('Color',[0.85, 0.33, 0.1],'Marker','s','Name','FRU'), ...
    'galactose', struct('Color',[0.93, 0.69, 0.13],'Marker','d','Name','GAL'), ...
    'xylose',    struct('Color',[0.49, 0.18, 0.56],'Marker','^','Name','XYL'), ...
    'lactose',   struct('Color',[0.47, 0.67, 0.19],'Marker','v','Name','LAC'), ...
    'sucrose',   struct('Color',[0.64, 0.08, 0.18],'Marker','p','Name','SUC') ...
);

%% Sce
% Load table and remove rows with missing values in experimental/predicted rates
T = readtable('CSource_res.xlsx', 'Sheet', 'Sce');
rows_with_content = ~ismissing(T{:,6});
T(rows_with_content, :) = [];
carbon_sources = T.("Carbon_source");
mu_exp = T.("mu_in_vivo_h_1_");
mu_pred = T.("mu_in_silico_h_1_");

% Define plotting styles (color, marker, label) for carbon sources
styles = struct( ...
    'glucose',   struct('Color',[0, 0.45, 0.74],'Marker','o','Name','GLC'), ...
    'sucrose',   struct('Color',[0.64, 0.08, 0.18],'Marker','p','Name','SUC'), ...
    'maltose',   struct('Color',[1.0, 0.0, 0.5],'Marker','<','Name','MAL'), ...
    'galactose', struct('Color',[0.93, 0.69, 0.13],'Marker','d','Name','GAL'), ...
    'fructose',  struct('Color',[0.85, 0.33, 0.1],'Marker','s','Name','FRU') ...
);

%% Ppa
% Load table and remove rows with missing values in experimental/predicted rates
T = readtable('CSource_res.xlsx', 'Sheet', 'Ppa');
rows_with_content = ~ismissing(T{:,6});
T(rows_with_content, :) = [];
carbon_sources = T.("Carbon_source");
mu_exp = T.("mu_in_vivo_h_1_");
mu_pred = T.("mu_in_silico_h_1_");

% Define plotting styles (color, marker, label) for carbon sources
styles = struct( ...
    'glucose',   struct('Color',[0, 0.45, 0.74],'Marker','o','Name','GLC'), ...
    'glycerol',  struct('Color',[1.0, 0.5, 0.0],'Marker','h','Name','GLY'), ...
    'methanol',  struct('Color',[0.0, 0.75, 0.75],'Marker','>','Name','MEOH') ...
);


%%  Plotting
figure;
hold on;

set(gca, 'FontName', 'Arial', 'FontSize', 7,'FontWeight', 'bold');
% Plot experimental vs. predicted growth rates for each carbon source
unique_sources = unique(carbon_sources);
for i = 1:length(unique_sources)
    src = lower(unique_sources{i});  
    idx = strcmpi(carbon_sources, unique_sources{i});

    if isfield(styles, src)
        h = scatter(mu_exp(idx), mu_pred(idx), 20, styles.(src).Color, 'filled', styles.(src).Marker, 'DisplayName', styles.(src).Name);

        legend_handles(end+1) = h;
        legend_names{end+1} = styles.(src).Name;
    else
        h = scatter(mu_exp(idx), mu_pred(idx), 20,'k', 'filled', 'o', 'DisplayName', unique_sources{i});

        legend_handles(end+1) = h;
        legend_names{end+1} = unique_sources{i};
    end
end


% Compute correlation metrics
all_exp = mu_exp;
all_pred = mu_pred;

% Pearson correlation coefficient and p-value
[pearson_r, p_val] = corr(all_exp, all_pred, 'Type', 'Pearson');
p = polyfit(all_exp, all_pred, 1);  
fit_line = polyval(p, all_exp);     
% Reference line (y = x)
max_val = max([all_exp;all_pred]) * 1.2;  
plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

% Axis labels and formatting
xlim([0 max_val]);
ylim([0 max_val]);
xlabel('Experimental growth rate [gDW/h]', 'FontSize', 7, 'FontWeight', 'bold');
ylabel('Predicted growth rate [gDW/h]', 'FontSize', 7, 'FontWeight', 'bold');
set(gcf,'units', 'centimeters','position',[10 10 5 5],'Color','none');
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1 1 3.5 3.5],'Color','none');
% Annotation with r and p-value
text(0.1 * max_val, 0.9 * max_val, sprintf('r = %.2f', pearson_r), 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
exp_val = floor(log10(p_val));
mantissa = p_val / 10^exp_val;
exp_str = num2str(abs(exp_val));
sup_map = '⁰¹²³⁴⁵⁶⁷⁸⁹';  
sup_digits = sup_map(double(exp_str) - double('0') + 1);
if exp_val < 0
    sup_digits = ['⁻' sup_digits];
else
    sup_digits = ['⁺' sup_digits];
end

p_text = sprintf('p = %.2f × 10%s', mantissa, sup_digits);
text(0.1 * max_val, 0.8 * max_val, p_text, 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
text(0.1 * max_val, 0.7 * max_val, 'Ppa', 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');

[unique_names, ia] = unique(legend_names, 'stable');  
unique_handles = legend_handles(ia);

legend(unique_handles, unique_names, 'Location', 'bestoutside', 'FontSize', 7, 'Box', 'off');
hold off;