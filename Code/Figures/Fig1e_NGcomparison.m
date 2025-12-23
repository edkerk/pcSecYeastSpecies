%% Fig1e_NGcomparison
% Read NG number information from three yeast models
% Each dataset corresponds to one strain: P. pastoris (Ppa), K. marxianus (Kma), S. cerevisiae (Sce)
addpath('../../Data/')
% K. marxianus
[~,~,protein_infoKMX] = xlsread('protein_information_KM.xlsx');
NG_KMX = protein_infoKMX(2:end,6);

% S. cerevisiae
[~,~,protein_infoSCE] = xlsread('Protein_Information_SCE.xlsx');
NG_SCE = protein_infoSCE(2:end,6);

% P. pastoris
[~,~,protein_infoPP] = xlsread('Protein_Information_PP.xlsx');
NG_PP = protein_infoPP(2:end,6);

%% Plot Violin Plots
% Convert cell to numeric, filter out empty entries and nonzero values
clean2num      = @(c) cell2mat(c(~cellfun(@isempty, c)));
filter_nonzero = @(x) x(x > 0);

data.Ppa = filter_nonzero(clean2num(NG_PP));
data.Kma = filter_nonzero(clean2num(NG_KMX));
data.Sce = filter_nonzero(clean2num(NG_SCE));
colors = [141 188 208; 110 176 149; 227 132 147]./255;  

figure; hold on;
set(gca, 'Color', 'none', 'FontName', 'Arial', ...
    'FontSize', 7, 'FontWeight', 'bold');

fields = fieldnames(data);

for i = 1:length(fields)
    % Extract strain data
    strain_data = data.(fields{i});
    strain_data = strain_data(~isnan(strain_data));

    % Define main and derived colors
    main_color    = colors(i,:);
    box_color     = min(main_color + 0.15, 1);
    median_color  = max(main_color - 0.15, 0);
    outlier_color = min(main_color + 0.25, 1); 
    
    % Kernel density estimation for violin shape
    [f, yi] = ksdensity(strain_data);
    f = f * 3;  % adjust width scaling

    % Draw violin
    fill([f, -f(end:-1:1)]*0.5 + i, [yi, yi(end:-1:1)], ...
         main_color, 'FaceAlpha', 0.4, ...
         'EdgeColor', main_color, 'LineWidth', 1);

    % Outliers (by quartile method)
    outli = isoutlier(strain_data, 'quartiles');
    scatter(i*ones(sum(outli),1), strain_data(outli), 10, ...
            'MarkerFaceColor', outlier_color, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

    % Quartiles and median
    qt25 = quantile(strain_data, 0.25);
    qt75 = quantile(strain_data, 0.75);
    med  = median(strain_data);

    % Whiskers
    plot([i i], [min(strain_data(~outli)) qt25], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    plot([i i], [qt75 max(strain_data(~outli))], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);

    % Box 
    fill(i + 0.1*[-1 1 1 -1], [qt25 qt25 qt75 qt75], ...
         box_color, 'FaceAlpha', 0.9, ...
         'EdgeColor', main_color, 'LineWidth', 1.2);

    % Median line
    plot(i + 0.1*[-1 1], [med med], 'Color', median_color, 'LineWidth', 2.5);
end

%% Axis Formatting
set(gca,'units', 'centimeters', ...
    'XTick', 1:3, 'XTickLabel', {'Ppa','Kma','Sce'}, ...
    'LineWidth', 0.5, 'Position', [1 0.5 3 4]);

set(gcf, 'units', 'centimeters', 'position', [10 10 4.5 5]); 

ylabel('NG Number', 'FontSize', 7);

ylim([0 max(structfun(@max, data)) * 1.2]);
xlim([0.5 3.5]);
