%% Fig3f_EnzymeTempSensitivity
% Visualize temperature sensitivity of ETC enzymes across yeasts
addpath('../../Results/')
filePath = 'TDM_pathway.xlsx';
opts = detectImportOptions(filePath);
data = readtable(filePath, opts);

pathways = data.pathway;
values = [data.KMX, data.SCE, data.PP];


%% Plot heatmap
figure;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

maxValue = 12;
customMap = [linspace(1,1,256)' linspace(1,0,256)' linspace(1,0,256)'];

h = heatmap(pathways, {'Kma', 'Sce', 'Ppa'}, values', ...
    'Colormap', customMap, ...
    'ColorLimits', [0 maxValue], ...
    'ColorbarVisible', 'on', ...
    'CellLabelColor', 'none');

h.XLabel = 'Pathways';
h.YLabel = 'Species';
h.GridVisible = 'off'; 

% Formatting
s = struct(h); 
ax = s.Axes;
ax.XTickLabelRotation = 30;
ax.XAxis.FontSize = 7;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

% Figure size and axis position
set(gcf, 'units', 'centimeters', 'position', [10 10 16 6.6]); 
set(gca, 'units', 'centimeters', 'Position', [2.5 4.5 12 2]);

