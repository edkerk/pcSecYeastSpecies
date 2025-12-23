%% FigS2abc_TPgrowth
addpath('../../Results/Growth_rate_TP/')

% Define target proteins and short labels
targetProteins = {'galactosidase';'Glucoseoxidase';'Inulinase';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'Humantransferin';'Hemoglobin';'PHO';'Amylase';'Insulin';'BGL'};
 targetProteins_rb = {'Gal';'GOX';'Inu';'SOD';'EST';'HSA';'COL';'PH20';'hLF';'TRF';'Hb';'PHO';'AMS';'IP';'BGL'};
color = [31,119,180; 214,39,40; 44,160,44; 148,103,189; 255,127,14; 227,119,194; 127,127,127; 188,189,34; 23,190,207; 140,86,75; 255,152,150; 174,199,232; 152,223,138; 196,156,148; 247,182,210]./255;

%% KMX

figure
hold on
set(gca, 'FontName', 'Arial', 'FontSize', 7,'FontWeight', 'bold');

mu_max = 0;
q_tp_max = 0;
% Loop over target proteins and plot μ–q curves
TP = targetProteins_rb;
for i = 1:length(targetProteins)
    display([num2str(i) '/' num2str(length(targetProteins))]);
    load(['modelKMX_',targetProteins{i},'.mat']);
    load(['fluxesKMX_',targetProteins{i},'.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    % Extract growth rate (μ) and target protein exchange flux
    tp_ex = [targetProteins{i},' exchange'];
    mu = round(fluxes(ismember(model.rxns,'r_1913'),:),2);
    [mu,b] = sort(mu,'ascend');
    fluxes = fluxes(:,b);
    q_tp = fluxes(ismember(model.rxns,tp_ex),:);
    mu_max   = max(mu_max, max(mu));
    q_tp_max = max(q_tp_max, max(q_tp));
    % Plot production vs growth rate
    plot(mu,q_tp,'-','LineWidth',1,'Color',color(i,:))

end

% Figure formatting
set(gcf,'units', 'centimeters','position',[10 10 5 5],'Color','none');
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1.5 1 3 3.5],'XLim', [0 mu_max*1.1], 'YLim', [0 q_tp_max*1.1],'Color','none');

yticks_original = get(gca, 'YTick');
scale_exp = -3;  
yticks_scaled = yticks_original / 10^scale_exp;

yticklabels(arrayfun(@(x) sprintf('%.0f', x), yticks_scaled, 'UniformOutput', false));
% Add scale annotation
text_pos_x = -0.1 * mu_max;  
text_pos_y = q_tp_max * 1.2; 
text(text_pos_x, text_pos_y, ['×10^{', num2str(scale_exp), '}'], ...
    'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
text(0.8*mu_max, 0.9*q_tp_max, 'Kma',...
    'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
% Axis labels
% % legend(string(TP),'Fontsize',6,'Box','off')
ylabel({'Protein production rate'; '[mmol/gDW/h]'},'FontSize',7);
xlabel('Growth rate [1/h]','FontSize',7);


%% SCE

figure
hold on
set(gca, 'FontName', 'Arial', 'FontSize', 7,'FontWeight', 'bold');

mu_max = 0;
q_tp_max = 0;
% Loop over target proteins and plot μ–q curves
TP = targetProteins;
for i = 1:length(targetProteins)
    display([num2str(i) '/' num2str(length(targetProteins))]);
    load(['modelSCE_',targetProteins{i},'.mat']);
    load(['fluxesSCE_',targetProteins{i},'.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    % Extract growth rate (μ) and target protein exchange flux
    tp_ex = [targetProteins{i},' exchange'];
    mu = round(fluxes(ismember(model.rxns,'r_2111'),:),2);
    [mu,b] = sort(mu,'ascend');
    fluxes = fluxes(:,b);
    q_tp = fluxes(ismember(model.rxns,tp_ex),:);
    mu_max   = max(mu_max, max(mu));
    q_tp_max = max(q_tp_max, max(q_tp));
    % Plot production vs growth rate
    plot(mu,q_tp,'-','LineWidth',1,'Color',color(i,:))

end

% Figure formatting
set(gcf,'units', 'centimeters','position',[10 10 5 5],'Color','none');
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1.5 1 3 3.5],'XLim', [0 mu_max*1.1], 'YLim', [0 q_tp_max*1.1],'Color','none');

yticks_original = get(gca, 'YTick');
scale_exp = -3;  
yticks_scaled = yticks_original / 10^scale_exp;

yticklabels(arrayfun(@(x) sprintf('%.0f', x), yticks_scaled, 'UniformOutput', false));
% Add scale annotation
text_pos_x = -0.1 * mu_max;  
text_pos_y = q_tp_max * 1.2; 
text(text_pos_x, text_pos_y, ['×10^{', num2str(scale_exp), '}'], ...
    'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
text(0.8*mu_max, 0.9*q_tp_max, 'Sce',...
    'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
% Axis labels
ylabel({'Protein production rate'; '[mmol/gDW/h]'},'FontSize',7);
xlabel('Growth rate [1/h]','FontSize',7);



%% PP

figure
hold on
set(gca, 'FontName', 'Arial', 'FontSize', 7,'FontWeight', 'bold');

mu_max = 0;
q_tp_max = 0;
TP = extractAfter(targetProteins, "_");
TP = extractAfter(TP, "_");

for i = 1:length(targetProteins)
    display([num2str(i) '/' num2str(length(targetProteins))]);
    load(['modelPP_',targetProteins{i},'.mat']);
    load(['fluxesPP_',targetProteins{i},'.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    % Extract growth rate (μ) and target protein exchange flux
    tp_ex = [targetProteins{i},' exchange'];
    mu = round(fluxes(ismember(model.rxns,'BIOMASS'),:),4);
    [mu,b] = sort(mu,'ascend');
    fluxes = fluxes(:,b);
    q_tp = fluxes(ismember(model.rxns,tp_ex),:);
    mu_max   = max(mu_max, max(mu));
    q_tp_max = max(q_tp_max, max(q_tp));
    % Plot production vs growth rate
    plot(mu,q_tp,'-','LineWidth',1,'Color',color(i,:))

end

% Figure formatting
set(gcf,'units', 'centimeters','position',[10 10 5 5],'Color','none');
set(gca,'units', 'centimeters', 'LineWidth', 0.5,'Position', [1.5 1 3 3.5],'XLim', [0 mu_max*1.1], 'YLim', [0 q_tp_max*1.1],'Color','none');

yticks_original = get(gca, 'YTick');
scale_exp = -3;  
yticks_scaled = yticks_original / 10^scale_exp;

yticklabels(arrayfun(@(x) sprintf('%.0f', x), yticks_scaled, 'UniformOutput', false));
% Add scale annotation
text_pos_x = -0.1 * mu_max; 
text_pos_y = q_tp_max * 1.2; 
text(text_pos_x, text_pos_y, ['×10^{', num2str(scale_exp), '}'], ...
    'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'tex');
text(0.8*mu_max, 0.9*q_tp_max, 'Ppa',...
    'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
% Axis labels
ylabel({'Protein production rate'; '[mmol/gDW/h]'},'FontSize',7);
xlabel('Growth rate [1/h]','FontSize',7);