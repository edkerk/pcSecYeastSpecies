%% Fig1b_ModelComparisonPpa

% Data preparation
categories = {'Genes', 'Reactions', 'Metabolites'};
pcSecPchia = [1417, 29026, 20195];   % pcSecPchia model
iMT1026v3     = [1025, 2237, 1706];      % iMT1026v3 reference model

%% Figure setup
figure('Position', [100, 100, 700, 500]); 
set(groot, 'defaultAxesFontSize', 7);

% Define color scheme
colors = [207/255 218/255 220/255;   % iMT1026v3 
          141/255, 188/255, 208/255];  % pcSecPchia 

%% Bar plot construction
bar_positions = 1:3;
bar_width = 0.35;

% Plot iSM996 (left bars)
b1 = bar(bar_positions - bar_width/2, iMT1026v3, bar_width, ...
         'FaceColor', colors(1,:), 'EdgeColor','none');
hold on;

% Plot pcSecKmarx (right bars)
b2 = bar(bar_positions + bar_width/2, pcSecPchia, bar_width, ...
         'FaceColor', colors(2,:), 'EdgeColor','none');

%% Axis formatting
set(gca, 'XTick', bar_positions, 'XTickLabel', categories);
ylim([0 80000]);   % Full range
ytickformat('%.0f');

% Axis label
ylabel('Number', 'FontSize', 7);

%% Legend
legend({'iMT1026v3', 'pcSecPchia'}, ...
       'Box', 'off', 'FontSize', 7, ...
       'Location','northwest');  

%% Aesthetic settings
set(gca, 'FontSize', 7);
ax = gca;
ax.Box = 'off';
