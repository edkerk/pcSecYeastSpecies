%% FigS1cde_ESC
% Enzyme sensitivity analysis of K. marxianus, S. cerevisiae and K. phaffii across temperature gradients.
% Computes FCC matrix across temperature gradients and plots top 20 genes. 
%% Kma
addpath('../../Enzymedata/')
addpath('../../Model/')

% load model and param
load('enzymedata_KMX.mat')
load('enzymedataMachine_KMX.mat')
load('enzymedataSEC_KMX.mat')
load('enzymedataDummyER_KMX.mat')
load('pcSecKmarx.mat')


t_list = [20 25 30 35 40 45 50];
num_genes = length(enzymedata.proteins); 
num_temps = length(t_list);
fcc_matrix = zeros(num_genes, num_temps);

%  Read flux files
for k = 1:num_temps
    t = t_list(k);
    pattern = sprintf('fluxestuningkcat1.1dummyKM_%d_*_1_1382.mat', t);
    files = dir(pattern);
    if isempty(files)
        error('No flux file found for temperature %d°C', t);
    end
    filename = files(1).name;

    parts = split(filename, '_');
    mu_val = str2double(parts{3}) / 100;   
    load(filename);  
    fcc_matrix(:, k) = (sol_dummy - sol_dummy_ref) / (sol_dummy_ref * mu_val);
end

%% Sce
% load model and param
load('enzymedata_SCE.mat')
load('enzymedataMachine_SCE.mat')
load('enzymedataSEC_SCE.mat')
load('enzymedataDummyER_SCE.mat')
load('pcSecYeast.mat')

t_list = [20 25 30 35 40];
num_genes = length(enzymedata.proteins); 
num_temps = length(t_list);
fcc_matrix = zeros(num_genes, num_temps);

%  Read flux files
for k = 1:num_temps
    t = t_list(k);
    pattern = sprintf('fluxestuningkcat1.1dummySCE_%d_*_1_1382.mat', t);
    files = dir(pattern);
    if isempty(files)
        error('No flux file found for temperature %d°C', t);
    end
    filename = files(1).name;

    parts = split(filename, '_');
    mu_val = str2double(parts{3}) / 100;   
    load(filename);  
    fcc_matrix(:, k) = (sol_dummy - sol_dummy_ref) / (sol_dummy_ref * mu_val);
end


%% Ppa
% load model and param
load('enzymedata_PP.mat')
load('enzymedataMachine_PP.mat')
load('enzymedataSEC_PP.mat')
load('enzymedataDummyER_PP.mat')
load('pcSecPichia.mat')

t_list = [20 25 30 35 40];
num_genes = length(enzymedata.proteins); 
num_temps = length(t_list);
fcc_matrix = zeros(num_genes, num_temps);

%  Read flux files
for k = 1:num_temps
    t = t_list(k);
    pattern = sprintf('fluxestuningkcat1.1dummyPP_%d_*_1_1382.mat', t);
    files = dir(pattern);
    if isempty(files)
        error('No flux file found for temperature %d°C', t);
    end
    filename = files(1).name;

    parts = split(filename, '_');
    mu_val = str2double(parts{3}) / 100;   
    load(filename);  
    fcc_matrix(:, k) = (sol_dummy - sol_dummy_ref) / (sol_dummy_ref * mu_val);
end


%% Filter genes with zero FCC at highest temperature
highest_temp_idx = num_temps;
zero_genes = (fcc_matrix(:, highest_temp_idx) == 0);
filtered_fcc = fcc_matrix(~zero_genes, :);
filtered_genes = enzymedata.proteins(~zero_genes);

% Rank top 20 genes based on highest-T FCC
[~, sort_idx] = sort(filtered_fcc(:, highest_temp_idx), 'descend');
topN = min(20, length(filtered_genes));
top_genes = filtered_genes(sort_idx(1:topN));
top_fcc = filtered_fcc(sort_idx(1:topN), :);


% Colormap
n = 256;
redwhiteMap = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% Plot
figure;
set(gca, 'Color', 'none','FontName', 'Arial', 'FontSize', 7);

h = heatmap(t_list, top_genes, top_fcc,'GridVisible', 'off','FontName', 'Arial');
colormap(redwhiteMap);
caxis([0, 10]); 

h.XLabel = 'Temperature (°C)';
h.YLabel = 'Gene';
h.YDisplayLabels = strrep(top_genes, '_', '\_');

ax = struct(h).Axes;
ax.FontName = 'Arial';
ax.FontSize = 7;
% % ax.FontWeight = 'bold';
% % ax.XLabel.FontWeight = 'bold';
% % ax.YLabel.FontWeight = 'bold';


set(gcf, 'units', 'centimeters', 'position', [5 5 5 7.5],'Color', 'none');
set(gca, 'units', 'centimeters', 'Position', [1.5 1 2.5 6]);



