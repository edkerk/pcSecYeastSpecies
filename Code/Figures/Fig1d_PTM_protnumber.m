%% Fig1d_PTM_protnumber.m
% Compare PTM-related protein numbers across three yeasts
% Count proteins with non-zero entries in SP, TRMM, ThroughER, NG, OG, DSB, GPI, total protein count
% Species: P. pastoris (Ppa), K. marxianus (Kma), S. cerevisiae (Sce)
addpath('../../Data/')
%% Load data
% K. marxianus
[~,~,protein_infoKma] = xlsread('protein_information_KM.xlsx');
% S. cerevisiae
[~,~,protein_infoSce] = xlsread('Protein_Information_SCE.xlsx');
% P. pastoris
[~,~,protein_infoPpa] = xlsread('Protein_Information_PP.xlsx');

%% Define column indices (based on provided Excel structure)
col_ThroughER = 3;
col_SP        = 4;
col_DSB       = 5;
col_NG        = 6;
col_OG        = 7;
col_TRMM      = 8;
col_GPI       = 9;

labels_full = {'ThroughER','Signal peptide','Disulfide site','N-glycosylation',...
               'O-glycosylation','Transmembrane','GPI site','All protein'};
labels_short = {'ThroughER','SP','DSB','NG','OG','TRMM','GPI','All protein'}; 

%% Function: count non-zero entries per feature
count_features = @(data) [ ...
    sum(cellfun(@(x) isnumeric(x) && ~isnan(x) && x>0, data(2:end,col_ThroughER))), ...
    sum(cellfun(@(x) isnumeric(x) && ~isnan(x) && x>0, data(2:end,col_SP))), ...
    sum(cellfun(@(x) isnumeric(x) && ~isnan(x) && x>0, data(2:end,col_DSB))), ...
    sum(cellfun(@(x) isnumeric(x) && ~isnan(x) && x>0, data(2:end,col_NG))), ...
    sum(cellfun(@(x) isnumeric(x) && ~isnan(x) && x>0, data(2:end,col_OG))), ...
    sum(cellfun(@(x) isnumeric(x) && ~isnan(x) && x>0, data(2:end,col_TRMM))), ...
    sum(cellfun(@(x) isnumeric(x) && ~isnan(x) && x>0, data(2:end,col_GPI))), ...
    size(data,1)-1 ];  % total proteins

resKma = count_features(protein_infoKma);
resSce = count_features(protein_infoSce);
resPpa = count_features(protein_infoPpa);

res = [resPpa; resKma; resSce]; 
[~,order] = sort(sum(res,1),'ascend');
res_sorted = res(:,order);
labels_sorted = labels_short(order);

%% Plot bar chart
figure; hold on;
set(gcf,'Units','inches','Position',[1 1 8 5],'Color',[1 1 1]);

cols = [
    141,188,208;   % Ppa 
    110,176,149;   % Kma 
    227,132,147    % Sce 
]/255;

b = bar(res_sorted','grouped','BarWidth',0.7);
for i=1:3
    b(i).FaceColor = cols(i,:);
end

ax = gca;
ax.LineWidth = 0.6;
ax.FontSize = 12;
ax.FontName = 'Arial';
xticks(1:length(labels_sorted));
xticklabels(labels_sorted);
xtickangle(45);
ylabel('Protein number','FontSize',7,'FontWeight','bold');

legend({'pcSecPichia','pcSecKmarx','pcSecYeast'},...
    'FontSize',7,'FontName','Arial','Box','off','Location','northwest');

hold off;
