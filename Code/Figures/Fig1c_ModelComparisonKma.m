%% Fig1c_ModelComparisonKma

% Data preparation
categories = {'Genes', 'Reactions', 'Metabolites'};
pcSecKmarx = [1382, 49471, 28888];   % pcSecKmarx model
iSM996     = [996, 1913, 1531];      % iSM996 reference model

%% Figure setup
figure('Position', [100, 100, 700, 500]); 
set(groot, 'defaultAxesFontSize', 7);

% Define color scheme 
colors = [207/255 218/255 220/255;   % iSM996 
          110/255 176/255 149/255];  % pcSecKmarx 

%% Bar plot construction
bar_positions = 1:3;
bar_width = 0.35;

% Plot iSM996 (left bars)
b1 = bar(bar_positions - bar_width/2, iSM996, bar_width, ...
         'FaceColor', colors(1,:), 'EdgeColor','none');
hold on;

% Plot pcSecKmarx (right bars)
b2 = bar(bar_positions + bar_width/2, pcSecKmarx, bar_width, ...
         'FaceColor', colors(2,:), 'EdgeColor','none');

%% Axis formatting
set(gca, 'XTick', bar_positions, 'XTickLabel', categories);
ylim([0 80000]);   % Full range
ytickformat('%.0f');

% Axis label
ylabel('Number', 'FontSize', 7);

%% Legend
legend({'iSM996', 'pcSecKmarx'}, ...
       'Box', 'off', 'FontSize', 7, ...
       'Location','northwest');  

%% Aesthetic settings
set(gca, 'FontSize', 7);
ax = gca;
ax.Box = 'off';
