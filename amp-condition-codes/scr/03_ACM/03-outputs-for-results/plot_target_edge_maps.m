% Visualizes how edges/networks are targeted and modulated by the rest of the brain.
% Supports both 'Network' (17x17) and 'Node' (200x200) resolution.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST2_LR'; 
group_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Group_Level', session);

% --- CHOOSE RESOLUTION LEVEL HERE ---
% Set to 'Network' for 17x17 atlas mapping, or 'Node' for 200x200 raw mapping
display_level = 'Network'; % Options: 'Network', 'Node'
t_threshold = 2.5; % Significance threshold for T-Values

% Load Atlas Info (Required for both to know matrix sizes)
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17;
n_rois = 200;

% Edge mapping indices
[u_map, v_map] = find(triu(true(n_rois), 1));
n_edges = length(u_map);

%% 2. Load Group-Level Matrices [19900 x 200]
fprintf('Loading Group-Level matrices...\n');
V = load(fullfile(group_dir, 'Group_Variance_Consensus.mat')); var_mat = V.Consensus_Percent;
T = load(fullfile(group_dir, 'Group_Slope_Tmap.mat'));         t_mat = T.T_map;
Z = load(fullfile(group_dir, 'Group_ZCR_Median.mat'));         zcr_mat = Z.Median;
P = load(fullfile(group_dir, 'Group_Poly_Prevalence.mat'));    poly_mat = P.Mode;

%% 3. Data Aggregation (Target Edge Perspective: Reduce 200 modulators to 1 value per edge)
fprintf('Calculating Metrics per Edge...\n');

% --- A. Variance Metrics ---
var_mean = mean(var_mat, 2, 'omitnan');

% Entropy: Normalize to probability, H = -sum(p * log2(p))
p_var = var_mat ./ sum(var_mat, 2, 'omitnan');
p_var(p_var == 0 | isnan(p_var)) = 1; % log2(1) = 0, cancels out safely
var_entropy = -sum(p_var .* log2(p_var), 2);

% Top 5 Max Mean
var_sorted = sort(var_mat, 2, 'descend', 'MissingPlacement', 'last');
var_top5 = mean(var_sorted(:, 1:5), 2, 'omitnan');

% --- B. Slope (T-Value) Metrics ---
slope_pos_count = sum(t_mat > t_threshold, 2);
slope_neg_count = sum(t_mat < -t_threshold, 2);
slope_abs_mean  = mean(abs(t_mat), 2, 'omitnan');

% --- C. ZCR/MCR Metric ---
zcr_mean = mean(zcr_mat, 2, 'omitnan');

% --- D. Polynomial Degree Metrics (Distribution of Modulators) ---
poly_static  = (sum(poly_mat == 0, 2) ./ n_rois) * 100; % % of modulators causing Static
poly_linear  = (sum(poly_mat == 1, 2) ./ n_rois) * 100; % % causing Linear
poly_nonlin  = (sum(ismember(poly_mat, [2, 3]), 2) ./ n_rois) * 100; % % causing Non-linear

%% 4. Spatial Mapping Function
% Helper inline function to map [19900x1] to either [17x17] or [200x200]
map_to_space = @(edge_vec) map_edges(edge_vec, display_level, u_map, v_map, roi_nets, unique_nets, n_rois, n_nets);

% Map all vectors
map_v_mean = map_to_space(var_mean);
map_v_ent  = map_to_space(var_entropy);
map_v_top5 = map_to_space(var_top5);

map_s_pos  = map_to_space(slope_pos_count);
map_s_neg  = map_to_space(slope_neg_count);
map_s_abs  = map_to_space(slope_abs_mean);

map_z_mean = map_to_space(zcr_mean);

map_p_stat = map_to_space(poly_static);
map_p_lin  = map_to_space(poly_linear);
map_p_non  = map_to_space(poly_nonlin);

%% 5. Visualization
labels = unique_nets;
if strcmp(display_level, 'Node'), labels = []; end % Hide labels for Node level

% --- Figure 1: Variance Metrics ---
figure('Name', 'Variance Modulatory Impact', 'Position', [50, 50, 1500, 450], 'Color', 'w');
sgtitle(['Variance Sensitivity (' display_level ' Level)'], 'FontSize', 16, 'FontWeight', 'bold');
plot_subplot(1, 2, 1, map_v_mean, 'Mean Variance Susceptibility', labels);
plot_subplot(1, 2, 2, map_v_ent,  'Entropy of Modulators (Diversity)', labels);
% plot_subplot(1, 3, 3, map_v_top5, 'Mean of Top-5 Modulators (Specific Impact)', labels);

% --- Figure 2: Slope (T-Value) Metrics ---
figure('Name', 'Slope (Directional) Impact', 'Position', [100, 100, 1500, 450], 'Color', 'w');
sgtitle(['Slope / T-Value Sensitivity (' display_level ' Level)'], 'FontSize', 16, 'FontWeight', 'bold');
plot_subplot(1, 2, 1, map_s_pos, 'Count of Significant Facilitators (T > 2.5)', labels);
plot_subplot(1, 2, 2, map_s_neg, 'Count of Significant Inhibitors (T < -2.5)', labels);
% plot_subplot(1, 3, 3, map_s_abs, 'Mean Absolute T-Value', labels);

% --- Figure 3: ZCR/MCR Metric ---
figure('Name', 'ZCR/MCR Impact', 'Position', [150, 150, 600, 500], 'Color', 'w');
sgtitle(['ZCR / MCR Volatility (' display_level ' Level)'], 'FontSize', 16, 'FontWeight', 'bold');
plot_subplot(1, 1, 1, map_z_mean, 'Mean Zero/Mean-Crossing Rate', labels);

% --- Figure 4: Polynomial Degree Metrics ---
figure('Name', 'Polynomial Degree Profile', 'Position', [200, 200, 1500, 450], 'Color', 'w');
sgtitle(['Polynomial Degree Tendency (' display_level ' Level)'], 'FontSize', 16, 'FontWeight', 'bold');
plot_subplot(1, 3, 1, map_p_stat, 'Static Tendency (% Deg 0)', labels);
plot_subplot(1, 3, 2, map_p_lin,  'Linear Tendency (% Deg 1)', labels);
plot_subplot(1, 3, 3, map_p_non,  'Non-Linear Tendency (% Deg 2 & 3)', labels);


%% --- HELPER FUNCTIONS ---

function out_mat = map_edges(edge_vec, level, u, v, roi_nets, unique_nets, n_nodes, n_nets)
    if strcmp(level, 'Node')
        % Map directly to 200x200 symmetric matrix
        out_mat = zeros(n_nodes, n_nodes);
        idx = sub2ind([n_nodes, n_nodes], u, v);
        out_mat(idx) = edge_vec;
        out_mat = out_mat + out_mat'; % Make symmetric
        % Fill diagonal with NaN for better visualization
        out_mat(1:n_nodes+1:end) = NaN; 
    else
        % Map to 17x17 Network matrix
        sum_scores = zeros(n_nets, n_nets);
        count_scores = zeros(n_nets, n_nets);
        
        net_idx_map = zeros(n_nodes, 1);
        for i = 1:n_nodes
            net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i}));
        end
        
        for e = 1:length(u)
            ni = u(e); nj = v(e);
            net_i = net_idx_map(ni);
            net_j = net_idx_map(nj);
            val = edge_vec(e);
            
            if ~isnan(val)
                sum_scores(net_i, net_j) = sum_scores(net_i, net_j) + val;
                count_scores(net_i, net_j) = count_scores(net_i, net_j) + 1;
                
                if net_i ~= net_j
                    sum_scores(net_j, net_i) = sum_scores(net_j, net_i) + val;
                    count_scores(net_j, net_i) = count_scores(net_j, net_i) + 1;
                end
            end
        end
        out_mat = sum_scores ./ count_scores;
    end
end

function plot_subplot(r, c, idx, mat, tit, labels)
    subplot(r, c, idx);
    imagesc(mat, 'AlphaData', ~isnan(mat)); % Handle NaNs
    colormap(jet); colorbar;
    
    % Axis Setup
    axis square;
    title(tit, 'FontSize', 12);
    
    if ~isempty(labels)
        % Network Level settings
        xticks(1:length(labels)); yticks(1:length(labels));
        xticklabels(labels); yticklabels(labels);
        xtickangle(45);
        set(gca, 'TickDir', 'out', 'FontSize', 8);
        
        % Overlay Values
        max_v = max(mat(:), [], 'omitnan');
        min_v = min(mat(:), [], 'omitnan');
        for i = 1:size(mat, 1)
            for j = 1:size(mat, 2)
                val = mat(i, j);
                if ~isnan(val)
                    % Text color logic
                    norm_v = (val - min_v) / (max_v - min_v + eps);
                    if norm_v < 0.25 || norm_v > 0.8, col='w'; else, col='k'; end
                    
                    % Format text based on magnitude
                    if val > 10, txt = sprintf('%.1f', val); else, txt = sprintf('%.2f', val); end
                    text(j, i, txt, 'HorizontalAlignment', 'center', 'Color', col, 'FontSize', 8, 'FontWeight', 'bold');
                end
            end
        end
    else
        % Node Level Settings
        set(gca, 'XTick', [], 'YTick', []); % Hide ticks for 200x200
    end
end