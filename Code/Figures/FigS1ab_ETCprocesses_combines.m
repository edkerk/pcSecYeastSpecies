
%% FigS1ab_ETCprocesses_combines
addpath('../../Results/Temperature-sensitive_parameters_analysis/')


t_list = [10 15 20 25 30 35 40 45];

%% Sce
% Extract growth rates
% Maximum μ [h^-1] across temperatures for each constraint model
mu_lists = struct();
load("fluxes_SCE_allprocesses.mat")
model_names = fieldnames(Results);
for m = 1:length(model_names)
    model_name = model_names{m};
    data_store = Results.(model_name);
    
    mu_list = nan(1, length(t_list));
    
    for t = 1:length(t_list)
        temp_val = t_list(t);
        rows = data_store(data_store(:,1) == temp_val, :);
        rows = rows(~isnan(rows(:,2)) & all(~isnan(rows(:,3:5)),2), :);
        if ~isempty(rows) 
            mu_list(t) = max(rows(:,2));
        else
            mu_list(t) = 0; 
        end
    end
    
    mu_lists.(model_name) = mu_list;
end
% Relative μ normalized to reference maximum
disp('t_list =');
disp(t_list);
for m = 1:length(model_names)
    mu_lists_norm.(model_names{m}) = mu_lists.(model_names{m})/max(mu_lists.PCSEC);
    fprintf('mu_list_%s = %s;\n', model_names{m}, mat2str(mu_lists.(model_names{m}), 3));
    fprintf('mu_list_norm_%s = %s;\n', model_names{m}, mat2str(mu_lists_norm.(model_names{m}), 3));
end

%% Ppa
% Extract growth rates
% Maximum μ [h^-1] across temperatures for each constraint model
mu_lists = struct();
load("fluxes_PP_allprocesses.mat")
model_names = fieldnames(Results);
for m = 1:length(model_names)
    model_name = model_names{m};
    data_store = Results.(model_name);
    
    mu_list = nan(1, length(t_list));
    
    for t = 1:length(t_list)
        temp_val = t_list(t);
        rows = data_store(data_store(:,1) == temp_val, :);
        rows = rows(~isnan(rows(:,2)) & all(~isnan(rows(:,3:5)),2), :);
        if ~isempty(rows) 
            mu_list(t) = max(rows(:,2));
        else
            mu_list(t) = 0; 
        end
    end
    
    mu_lists.(model_name) = mu_list;
end
% Relative μ normalized to reference maximum
disp('t_list =');
disp(t_list);
for m = 1:length(model_names)
    mu_lists_norm.(model_names{m}) = mu_lists.(model_names{m})/max(mu_lists.PCSEC);
    fprintf('mu_list_%s = %s;\n', model_names{m}, mat2str(mu_lists.(model_names{m}), 3));
    fprintf('mu_list_norm_%s = %s;\n', model_names{m}, mat2str(mu_lists_norm.(model_names{m}), 3));
end

%% Plot
color = [0.75,0.35,0.40;0.85,0.65,0.45;0.55,0.70,0.60;0.45,0.65,0.72;0.58,0.48,0.65];

figure;
hold on;  
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold');

plot(t_list, mu_lists_norm.PCSEC, '-s', 'LineWidth', 2, 'MarkerSize', 3, 'Color', color(1,:));
plot(t_list, mu_lists_norm.NGAM, '-s', 'LineWidth', 2, 'MarkerSize', 3, 'Color', color(2,:));
plot(t_list, mu_lists_norm.kcatT, '-s', 'LineWidth', 2, 'MarkerSize', 3,  'Color', color(3,:));
plot(t_list, mu_lists_norm.kdegT, '-s', 'LineWidth', 2, 'MarkerSize', 3,  'Color', color(4,:));
plot(t_list, mu_lists_norm.ETC, '-s', 'LineWidth', 2, 'MarkerSize', 3,  'Color', color(5,:));
set(gcf, 'units', 'centimeters', 'position', [10 10 8 6]); 
set(gca, 'units', 'centimeters', 'LineWidth', 0.5, 'Position', [1 1 6 4]);


legend({'pcsec', 'pcsec+NGAM','pcsec+kcatT','pcsec+kdegT','pcsec+ETC'}, 'FontSize', 7, 'Location', 'northeast', 'Box', 'off');
ylabel('Specific Growth rate [1/h]', 'FontSize', 7);
xlabel('Temperature [°C]', 'FontSize', 7);

x_range = xlim;
y_range = ylim;

x_pos = x_range(1) + 0.05 * (x_range(2) - x_range(1));
y_pos = y_range(1) + 0.95 * (y_range(2) - y_range(1));

text(x_pos, y_pos, 'Sce', ...
    'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Margin', 1);

hold off;
