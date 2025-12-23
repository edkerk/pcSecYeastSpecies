%% Fig1f_SecSubsystems_protnumber.m
% Purpose: Compare client vs machinery proteins across subsystems
% Species: K. marxianus (Kma), S. cerevisiae (Sce), P. pastoris (Ppa)
% Data: model.mat + TableS1.xlsx (Secretory sheet)
% Output: Grouped stacked bar plot 
addpath('../../Data/')

%% Load data 
% Kma
load('pcSecKmarx.mat');  
model_kma = model;
[~,~,proteins_info_kmx] = xlsread('TableS1_KM.xlsx','Secretory');
proteins_info_kmx = proteins_info_kmx(2:end,[1,2,26]);

% Sce
load('pcSecYeast.mat');
model_sce = model;
[~,~,proteins_info_sce] = xlsread('TableS1_SCE.xlsx','Secretory');
proteins_info_sce = proteins_info_sce(2:end,[1,2,26]);

% Ppa
load('pcSecPichia.mat')
model_ppa = model;
[~,~,proteins_info_ppa] = xlsread('TableS1_PP.xlsx','Secretory');
proteins_info_ppa = proteins_info_ppa(2:end,[1,2,26]);

%% Identify subsystems
processes = unique(proteins_info_kmx(:,3));
resKma = zeros(length(processes),2);
resSce = zeros(length(processes),2);
resPpa = zeros(length(processes),2);

%% Count clients and machinery proteins
for i = 1:length(processes)
    % Kma
    resKma(i,1) = sum(strcmp(proteins_info_kmx(:,3), processes{i})); % Clients
    complex = proteins_info_kmx(strcmp(proteins_info_kmx(:,3), processes(i)),1);
    rxns = model_kma.rxns(endsWith(model_kma.rxns,complex));
    rxns = regexprep(rxns,'KLMA_','KLMA-'); 
    rxns = extractBefore(rxns,'_');
    rxns = regexprep(rxns,'KLMA-','KLMA_');
    resKma(i,2) = numel(unique(setdiff(rxns,{'dummyER','dummy'}))); % Machinery

    % Sce
    resSce(i,1) = sum(strcmp(proteins_info_sce(:,3), processes{i}));
    complex = proteins_info_sce(strcmp(proteins_info_sce(:,3), processes(i)),1);
    rxns = model_sce.rxns(endsWith(model_sce.rxns,complex));
    rxns = extractBefore(rxns,'_');
    resSce(i,2) = numel(unique(setdiff(rxns,{'dummyER','dummy'})));

    % Ppa
    resPpa(i,1) = sum(strcmp(proteins_info_ppa(:,3), processes{i}));
    complex = proteins_info_ppa(strcmp(proteins_info_ppa(:,3), processes(i)),1);
    rxns = model_ppa.rxns(endsWith(model_ppa.rxns,complex));
    for j = 1:length(rxns)
        [~,pos] = find(rxns{j}=='_',3);
        rxns{j} = extractBefore(rxns{j},pos(3));
    end
    resPpa(i,2) = numel(unique(rxns));
end

%% Sort subsystems by total proteins
[~,idx] = sort(sum(resKma+resSce+resPpa,2));
processes = processes(idx);
resKma = resKma(idx,:); 
resSce = resSce(idx,:); 
resPpa = resPpa(idx,:);

%% Plot: grouped stacked bar
figure;
hold on;
set(gcf,'Units','inches','Position',[1 1 8 5],'Color',[1 1 1]);

x = 1:length(processes); 
bw = 0.7; % bar width

% Colors
cols = {
    [ 99,132,146]/255, [141,188,208]/255;  % Ppa
    [ 64,114, 97]/255, [110,176,149]/255;  % Kma
    [161, 80, 92]/255, [227,132,147]/255   % Sce
};

% Draw bars
b1 = bar(x-0.2, resPpa, 'stacked','BarWidth',bw/3,'EdgeColor','k');
b2 = bar(x,     resKma, 'stacked','BarWidth',bw/3,'EdgeColor','k');
b3 = bar(x+0.2, resSce, 'stacked','BarWidth',bw/3,'EdgeColor','k');

for i=1:2
    b1(i).FaceColor = cols{1,i}; 
    b2(i).FaceColor = cols{2,i};
    b3(i).FaceColor = cols{3,i};
end

%% Axes and labels
ax = gca; 
ax.XLim = [0.5 length(processes)+0.5]; 
ax.YLim(1) = 0;
ax.TickDir = 'out'; 
ax.LineWidth = 0.6;
ax.FontSize = 7; 
ax.FontName = 'Arial';
xticks(x); 
xticklabels(processes); 
xtickangle(45);
ylabel('Protein number','FontSize',7,'FontWeight','bold');

%% Legend
legend({'Ppa Client','Ppa Machinery','Kma Client','Kma Machinery',...
        'Sce Client','Sce Machinery'}, ...
        'FontSize',7,'FontName','Arial','Box','off','Location','northwest');

hold off;
