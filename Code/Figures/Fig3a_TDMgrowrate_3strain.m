%% Fig3a_TDMgrowrate_3strain
% Compare the maximal growth rate of K. marxianus, S. cerevisiae,
% and K. phaffii across temperature gradients.
addpath('../../Results/Growth_rate_TDM/')

%% Kma
load('pcSecKmarx.mat')
t_list = [14:2:52];
mu_list = zeros(size(t_list)); 

for i = 1:length(t_list)
    display([num2str(i) '/' num2str(length(t_list))]);
    load(['fluxesKM_', num2str(t_list(i)), '.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    mu = round(fluxes(ismember(model.rxns,'r_1913'),:), 2);% biomass rxn
    
    if ~isempty(mu)
        max_mu = max(mu(mu > 0)); 
        mu_list(i) = max_mu;  
    else
        mu_list(i) = NaN;  
    end
end

%% Sce
load('pcSecYeast.mat')
t_list_ref = [10:2:42];

mu_list_ref = zeros(size(t_list_ref));
for i = 1:length(t_list_ref)
    display([num2str(i) '/' num2str(length(t_list_ref))]);
    load(['fluxesSCE_', num2str(t_list_ref(i)), '.mat']);
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    mu = round(fluxes(ismember(model.rxns,'r_2111'),:), 2);% biomass rxn
    if ~isempty(mu)
        max_mu = max(mu(mu > 0));  
        mu_list_ref(i) = max_mu;   
    else
        mu_list_ref(i) = NaN;  
    end
end

%% PP
load('pcSecPichia.mat')
t_list_PP = [10:2:44];

mu_list_PP = zeros(size(t_list_PP));
for i = 1:length(t_list_PP)
    display([num2str(i) '/' num2str(length(t_list_PP))]);
    load(['fluxesPP_', num2str(t_list_PP(i)), '.mat']); 
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);    
    mu = round(fluxes(ismember(model.rxns,'BIOMASS'),:), 2);% biomass rxn
    if ~isempty(mu)
        max_mu = max(mu(mu > 0));  
        mu_list_PP(i) = max_mu;   
    else
        mu_list_PP(i) = NaN;  
    end
end

%% Plot results
figure;
hold on;  
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');
% Plot temperature-dependent growth profiles
plot(t_list, mu_list, '-o', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', [0.7, 0.7, 1], 'Color', [110,176,149] / 255);
plot(t_list_ref, mu_list_ref, '-s', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', [1, 0.6, 0.6], 'Color', [227,132,147] / 255);
plot(t_list_PP, mu_list_PP, '-x', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', [1, 0.6, 0.6], 'Color', [141,188,208] / 255);
% Format figure
set(gcf, 'units', 'centimeters', 'position', [10 10 7 6]);
set(gca, 'units', 'centimeters', 'LineWidth', 0.5, 'Position', [1 1 5 4]);
ylim([0 1]);

legend({'Kma', 'Sce','Ppa'}, 'Location', 'northeast', 'FontSize', 7, 'Box', 'off');
ylabel('Max Growth rate [1/h]', 'FontSize', 7);
xlabel('Temperature [°C]', 'FontSize', 7);

hold off;
