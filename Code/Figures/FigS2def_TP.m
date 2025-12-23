%% FigS2deg_TP
%    This script categorizes and visualizes overexpression targets 
%    identified in three yeast species (K. marxianus, S. cerevisiae, 
%    and P. pastoris) according to their functional modules.
addpath('../../Results/FSEOF/')

%% All targets  
targetProteins_rb = {'Gal';'GOX';'Inu';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'TRF';'Hb';'PHO';'AMS';'IP';'BGL'};

%% Kma
% Load target gene lists
load("res_geneListKMX_M2_new.mat")
res_allgene = res_list(:);
res_allgene(cellfun(@isempty,res_allgene)) = [];
res_corted = unique(res_allgene);
res_corted(:,2) = repmat({'metabolic'},length(res_corted),1);
load('enzymedataSEC_KMX.mat')
secprotein = unique(enzymedataSEC.subunit);
% Reassign genes that belong to secretory components
[~,idx] = ismember(res_corted(:,1),secprotein);
res_corted(idx ~= 0,2) = {'Secretory'};
[~,~,proteins]=xlsread('TableS1_KM.xlsx','Secretory');
[~,idx] = ismember(res_corted(:,1),proteins(:,2));
res_corted(idx ~= 0,2) = proteins(idx(idx~=0),26);
res_corted(strcmp(res_corted(:,2),'CPY')|strcmp(res_corted(:,2),'ALP'),2) = {'Sorting'};
res_corted(strcmp(res_corted(:,2),'COPII')|strcmp(res_corted(:,2),'COPI'),2) = {'ER\_Golgi transport'};
cata_KMX = unique(res_corted(:,2));
% Count overexpression targets per category for each recombinant protein
cata_res_KMX = zeros(length(targetProteins), length(cata_KMX));
for i = 1:length(targetProteins)
    tmp = res_list(:,i);
    tmp(cellfun(@isempty,tmp)) = [];
    [~,idx] = ismember(tmp,res_corted(:,1));
    tmp_cata = tabulate(res_corted(idx,2));
    [~,idx2] = ismember(tmp_cata(:,1),cata_KMX);
    cata_res_KMX(i,idx2) = cell2mat(tmp_cata(:,2));
end


%% Sce
% Load target gene lists
load("res_geneListSCE_M2_new.mat")
res_allgene = res_list(:);
res_allgene(cellfun(@isempty,res_allgene)) = [];
res_corted = unique(res_allgene);
res_corted(:,2) = repmat({'metabolic'},length(res_corted),1);
load('enzymedataSEC_SCE.mat')
secprotein = unique(enzymedataSEC.subunit);
% Reassign genes that belong to secretory components
[~,idx] = ismember(res_corted(:,1),secprotein);
res_corted(idx ~= 0,2) = {'Secretory'};
[~,~,proteins]=xlsread('TableS1_SCE.xlsx','Secretory');
[~,idx] = ismember(res_corted(:,1),proteins(:,2));
res_corted(idx ~= 0,2) = proteins(idx(idx~=0),26);
res_corted(strcmp(res_corted(:,2),'CPY')|strcmp(res_corted(:,2),'ALP'),2) = {'Sorting'};
res_corted(strcmp(res_corted(:,2),'COPII')|strcmp(res_corted(:,2),'COPI'),2) = {'ER\_Golgi transport'};
cata_SCE = unique(res_corted(:,2));
% Count overexpression targets per category for each recombinant protein
cata_res_SCE = zeros(length(targetProteins), length(cata_SCE));
for i = 1:length(targetProteins)
    tmp = res_list(:,i);
    tmp(cellfun(@isempty,tmp)) = [];
    [~,idx] = ismember(tmp,res_corted(:,1));
    tmp_cata = tabulate(res_corted(idx,2));
    [~,idx2] = ismember(tmp_cata(:,1),cata_SCE);
    cata_res_SCE(i,idx2) = cell2mat(tmp_cata(:,2));
end


%% Ppa
% Load target gene lists
load("res_geneListPP_M2_new.mat")
targetProteins = cellfun(@(x) extractAfter(x, find(x == '_', 1, 'last')), targetProteins, 'UniformOutput', false);
res_allgene = res_list(:);
res_allgene(cellfun(@isempty,res_allgene)) = [];
res_corted = unique(res_allgene);
res_corted(:,2) = repmat({'metabolic'},length(res_corted),1);
load('enzymedataSEC_PP.mat')
secprotein = unique(enzymedataSEC.subunit);
% Reassign genes that belong to secretory components
[~,idx] = ismember(res_corted(:,1),secprotein);
res_corted(idx ~= 0,2) = {'Secretory'};
[~,~,proteins]=xlsread('TableS1_PP.xlsx','Secretory');
[~,idx] = ismember(res_corted(:,1),proteins(:,2));
res_corted(idx ~= 0,2) = proteins(idx(idx~=0),26);
res_corted(strcmp(res_corted(:,2),'CPY')|strcmp(res_corted(:,2),'ALP'),2) = {'Sorting'};
res_corted(strcmp(res_corted(:,2),'COPII')|strcmp(res_corted(:,2),'COPI'),2) = {'ER\_Golgi transport'};
cata_PP = unique(res_corted(:,2));
% Count overexpression targets per category for each recombinant protein
cata_res_PP = zeros(length(targetProteins), length(cata_PP));
for i = 1:length(targetProteins)
    tmp = res_list(:,i);
    tmp(cellfun(@isempty,tmp)) = [];
    [~,idx] = ismember(tmp,res_corted(:,1));
    tmp_cata = tabulate(res_corted(idx,2));
    [~,idx2] = ismember(tmp_cata(:,1),cata_PP);
    cata_res_PP(i,idx2) = cell2mat(tmp_cata(:,2));
end


%% Construct unified category set and color map
colorList = [0.70, 0.20, 0.85; 0.40, 0.30, 0.95; 0.10, 0.50, 0.95; 0.10, 0.70, 0.85; 0.10, 0.75, 0.30; 0.80, 0.85, 0.10; 1.00, 0.70, 0.10; 0.95, 0.45, 0.10; 0.85, 0.10, 0.10];
cata_all = unique([cata_KMX; cata_SCE; cata_PP]);
color_map = containers.Map(cata_all, num2cell(colorList(1:length(cata_all), :), 2));

 targetProteins_rb = {'Gal';'GOX';'Inu';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'TRF';'Hb';'PHO';'AMS';'IP';'BGL'};


%% Plot: K. marxianus
figure;
hold on
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
x = 1:length(targetProteins);
gap = 0.05;
x_new = x * gap;

cata_all_KMX = zeros(length(targetProteins), length(cata_all));
[~, idx_map_KMX] = ismember(cata_KMX, cata_all);
cata_all_KMX(:, idx_map_KMX) = cata_res_KMX;

b = bar(x_new, cata_all_KMX, 'stacked', 'LineWidth', 0.5, 'BarWidth', 0.4);
for i = 1:length(cata_all)
    b(i).FaceColor = color_map(cata_all{i});
    b(i).FaceAlpha = 0.6;
end

ylim([0,200]);
xlim([min(x_new)-0.05, max(x_new)+0.05]);
text(x_new(3), 180, 'Kma', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
ylabel('Overexpression target number','FontSize',7, 'FontWeight', 'bold');
xticks(x_new);
xticklabels(targetProteins_rb);
xtickangle(90);
set(gcf,'units', 'centimeters','position',[10 10 5 6]);
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1 2 3.8 3.5],'Color','none');

%% Plot: S. cerevisiae
figure;
hold on
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
x_new = x * gap;

cata_all_SCE = zeros(length(targetProteins), length(cata_all));
[~, idx_map_SCE] = ismember(cata_SCE, cata_all);
cata_all_SCE(:, idx_map_SCE) = cata_res_SCE;

b = bar(x_new, cata_all_SCE, 'stacked', 'LineWidth', 0.5, 'BarWidth', 0.4);
for i = 1:length(cata_all)
    b(i).FaceColor = color_map(cata_all{i});
    b(i).FaceAlpha = 0.6;
end

% % legend(cata_all, 'FontSize', 7, 'Location', 'northeastoutside');
ylim([0,200]);
xlim([min(x_new)-0.05, max(x_new)+0.05]);
text(x_new(3), 180, 'Sce', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
ylabel('Overexpression target number','FontSize',7, 'FontWeight', 'bold');
xticks(x_new);
xticklabels(targetProteins_rb);
xtickangle(90);
set(gcf,'units', 'centimeters','position',[10 10 5 6]);
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1 2 3.8 3.5],'Color','none');


%% Plot: P. pastoris
figure;
hold on
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
x_new = x * gap;

cata_all_PP = zeros(length(targetProteins), length(cata_all));
[~, idx_map_PP] = ismember(cata_PP, cata_all);
cata_all_PP(:, idx_map_PP) = cata_res_PP;

b = bar(x_new, cata_all_PP, 'stacked', 'LineWidth', 0.5, 'BarWidth', 0.4);000000000
for i = 1:length(cata_all)
    b(i).FaceColor = color_map(cata_all{i});
    b(i).FaceAlpha = 0.6;
end

lgd = legend(cata_all, 'FontSize', 7, 'Location', 'northeastoutside', 'Box', 'off');
lgd.ItemTokenSize = [20, 6];  
ylim([0,200]);
xlim([min(x_new)-0.05, max(x_new)+0.05]);
text(x_new(3), 180, 'Ppa', 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
ylabel('Overexpression target number','FontSize',7, 'FontWeight', 'bold');
xticks(x_new);
xticklabels(targetProteins_rb);
xtickangle(90);
set(gcf,'units', 'centimeters','position',[10 10 5 6]);
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1 2 3.8 3.5],'Color','none');