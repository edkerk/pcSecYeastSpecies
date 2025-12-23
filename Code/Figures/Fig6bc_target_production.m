%% Fig6bc_target_production
% Fold-change plotting of hLF average production rate (Welch t-test stars)
addpath('../../Results/Experimental_validation/')

% settings
fname      = 'hLF_overproduction_result.xlsx';
colNumber  = 'Number';
colGene    = 'Gene';                      
colYield   = {'foldchange1','foldchange2','foldchange3'};
colCat     = 'Metabolic/Secretory';
filterCategory = 'metabolic';             % 'secretory' | 'metabolic' | ''(no filter)

ieq = @(a,b) strcmpi(string(a), string(b));

%% Load and filter
T = readtable(fname,'PreserveVariableNames',true);
rm = false(height(T),1);
if ismember(colNumber, T.Properties.VariableNames), rm = rm | ieq(T.(colNumber),'hLF-BS'); end
if ismember(colGene,   T.Properties.VariableNames), rm = rm | ieq(T.(colGene)  ,'hLF-BS'); end
T(rm,:) = [];

% Category filter
if ~isempty(filterCategory)
    want   = ieq(T.(colCat), filterCategory);
    if strcmpi(filterCategory,'secretory')
        isCtrl = ieq(T.(colNumber), 'GEM000_1');
    else
        isCtrl = ieq(T.(colNumber), 'GEM000_2');
    end
    T = T(want | isCtrl, :);
end

% take replicates; keep rows with >=2 non-NaN replicates
Yall = T{:, colYield};
validRows = sum(~isnan(Yall),2) >= 2;
T = T(validRows,:);  Yall = Yall(validRows,:);

% control first (GEM000)
if strcmpi(filterCategory,'secretory')
    ctrlIdx = ieq(T.(colNumber),'GEM000_1');
else
    ctrlIdx = ieq(T.(colNumber),'GEM000_2');
end
ord = [find(ctrlIdx); find(~ctrlIdx)];
T = T(ord,:);  Yall = Yall(ord,:);

%% Mean, SD, and two-sample t-tests
mu = mean(Yall,2,'omitnan');
sd = std (Yall,0,2,'omitnan');
xc = Yall(1,~isnan(Yall(1,:))); 

pvals = nan(height(T),1);
vartype_used = strings(height(T),1);  

for i = 2:height(T)
    xi = Yall(i,~isnan(Yall(i,:)));
    if numel(xi) >= 2 && numel(xc) >= 2 
        % Test for equal variance
        [~,p_var] = vartest2(xi, xc);  
        if p_var < 0.05
            [~,pvals(i)] = ttest2(xi, xc, 'Vartype','unequal','Tail','both');
            vartype_used(i) = "unequal (Welch)";
        else
            [~,pvals(i)] = ttest2(xi, xc, 'Vartype','equal','Tail','both');
            vartype_used(i) = "equal (Student)";
        end
    end
end



% convert to fold change for plotting
mu_ctrl   = mu(1);
mu_fc     = mu / mu_ctrl;     
sd_fc     = sd / mu_ctrl;      
Yall_fc   = Yall / mu_ctrl;   

% Significance labels
stars = strings(height(T),1);
for i = 2:height(T)
    if ~isnan(pvals(i)) && mu(i) > mu_ctrl
        if pvals(i) < 0.001, stars(i) = "***";
        elseif pvals(i) < 0.01, stars(i) = "**";
        elseif pvals(i) < 0.05, stars(i) = "*";
        end
    end
end

%% Plotting 
labels = T.(colGene);
if ~iscell(labels), labels = cellstr(labels); end
lab = labels;
lab{1} = labels{1};           
for i = 2:numel(lab)
    lab{i} = ['{\it ' labels{i} '}'];   
end

x = 1:height(T);
figure;
ax = axes('FontName','Arial','FontSize',7,'FontWeight','normal',...
          'Box','off','LineWidth',0.5,'Color','none'); 
hold(ax,'on');
ax.TickLabelInterpreter = 'tex';
% Colors
colEdge = [0 0 0];
colCtrl = [0.70 0.70 0.70];
if strcmpi(filterCategory,'secretory')
    colFace = [186 113 110]./255;
else
    colFace = [135 179 219]./255;
end


b = bar(x, mu_fc, 0.72, 'FaceColor','flat','EdgeColor','none');
b.CData = repmat(colFace,height(T),1);  b.CData(1,:) = colCtrl;

errorbar(x, mu_fc, sd_fc, 'LineStyle','none','Color',colEdge,'LineWidth',0.5);

rng(1);
for i = 1:height(T)
    yi = Yall_fc(i,~isnan(Yall_fc(i,:)));
    scatter(i + (rand(size(yi))-0.5)*0.20, yi, 6, 'o', ...
        'MarkerEdgeColor',colEdge, 'MarkerFaceColor','w','LineWidth',0.5);
end

yTop = mu_fc + sd_fc;  yOffset = 0.05 * max(yTop);
for i = 1:height(T)
    if strlength(stars(i))>0
        text(x(i), yTop(i)+yOffset, stars(i), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'Color',colEdge,'FontSize',7);
    end
end

set(ax,'XTick',x,'XTickLabel',lab,'XTickLabelRotation',45,...
       'LineWidth',1,'Box','off','FontName','Arial','FontSize',7);
ylabel({'Fold change of  rhLF', 'average production rate'},'FontSize',7);
xlabel('Overexpression targets','FontSize',7);

if strcmpi(filterCategory,'secretory')
    % secretory
    ylim([0, 6]);
    set(gcf,'units','centimeters','position',[10 10 7.5 5.5],'Color','none');
    set(ax,'units','centimeters','LineWidth',0.5,'Position',[1.2 1.5 6 3.5], ...
        'Color','none','FontName','Arial','FontSize',7);
else
    % metabolic
    ylim([0, 6]);
    set(gcf,'units','centimeters','position',[10 10 7.5 5.5],'Color','none');
    set(ax,'units','centimeters','LineWidth',0.5,'Position',[1.2 1.5 6 3.5], ...
        'Color','none','FontName','Arial','FontSize',7);
end

