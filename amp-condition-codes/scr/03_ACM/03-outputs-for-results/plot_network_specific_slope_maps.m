% Generates 17 Target-Edge Heatmaps (17x17) based on the T-Value (Slope).
% Each figure shows how ONE specific Modulator Network facilitates/inhibits the rest of the brain.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR'; 

% Input / Output Paths
group_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Group_Level', session);
fig_out_dir = fullfile('02 - Figs', 'Slope_Network_Specific_Modulators');

% Create Output Directory if it doesn't exist
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

t_threshold = 2.5; % Significance threshold for T-Values

% Load Atlas Info
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17;
n_rois = 200;

% Edge mapping indices
[u_map, v_map] = find(triu(true(n_rois), 1));
n_edges = length(u_map);

% Network Index Map (Map 200 ROIs to 17 Networks)
net_idx_map = zeros(n_rois, 1);
for i = 1:n_rois
    net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i}));
end

%% 2. Load Group-Level Slope Data [19900 x 200]
fprintf('Loading Group-Level T-Map...\n');
T = load(fullfile(group_dir, 'Group_Slope_Tmap.mat')); 
t_mat = T.T_map;

%% 3. Iterate Over Each Modulator Network (1 to 17)
fprintf('Generating specific heatmaps for %d Modulator Networks...\n', n_nets);

for mod_net = 1:n_nets
    mod_net_name = unique_nets{mod_net};
    
    % Find which nodes (out of 200) belong to this specific modulator network
    mod_nodes = find(net_idx_map == mod_net);
    
    % Extract T-values ONLY for these modulators [19900 x N_nodes_in_network]
    t_subset = t_mat(:, mod_nodes);
    
    % --- A. Calculate Metrics for this Subset ---
    % Count how many nodes IN THIS NETWORK significantly facilitated/inhibited each edge
    slope_pos_count = sum(t_subset > t_threshold, 2);
    slope_neg_count = sum(t_subset < -t_threshold, 2);
    
    % --- B. Map to 17x17 Target Space ---
    map_s_pos = map_edges(slope_pos_count, u_map, v_map, net_idx_map, n_nets);
    map_s_neg = map_edges(slope_neg_count, u_map, v_map, net_idx_map, n_nets);
    
    % --- C. Visualization ---
    fig = figure('Name', sprintf('Modulator: %s', mod_net_name), ...
                 'Position', [100, 100, 1400, 600], 'Color', 'w', 'Visible', 'off'); % 'off' makes saving faster
             
    sgtitle(sprintf('Modulator Network: %s (Target Edge Sensitivity)', mod_net_name), ...
            'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
        
    % Subplot 1: Facilitators
    plot_subplot(1, 2, 1, map_s_pos, 'Average Count of Facilitators (T > 2.5)', unique_nets);
    
    % Subplot 2: Inhibitors
    plot_subplot(1, 2, 2, map_s_neg, 'Average Count of Inhibitors (T < -2.5)', unique_nets);
    
    % --- D. Save Figure ---
    % Clean the network name for file saving (remove spaces/special chars if any)
    safe_name = matlab.lang.makeValidName(mod_net_name); 
    save_name_png = fullfile(fig_out_dir, sprintf('Modulator_%02d_%s.png', mod_net, safe_name));
    
    % Export as High-Quality PNG
    print(fig, save_name_png, '-dpng', '-r300');
    
    % Close figure to save memory
    close(fig);
    
    fprintf('Saved: %s\n', save_name_png);
end

fprintf('\nSuccess! All 17 figures saved in:\n%s\n', fig_out_dir);


%% --- HELPER FUNCTIONS ---

function out_mat = map_edges(edge_vec, u, v, net_idx_map, n_nets)
    % Maps [19900x1] edge vector to [17x17] Network matrix
    sum_scores = zeros(n_nets, n_nets);
    count_scores = zeros(n_nets, n_nets);
    
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

function plot_subplot(r, c, idx, mat, tit, labels)
    subplot(r, c, idx);
    imagesc(mat, 'AlphaData', ~isnan(mat)); 
    colormap(jet); colorbar;
    
    axis square;
    title(tit, 'FontSize', 14);
    
    % Network Level settings
    xticks(1:length(labels)); yticks(1:length(labels));
    xticklabels(labels); yticklabels(labels);
    xtickangle(45);
    
    % FIX: Use 'TickLabelInterpreter' instead of 'Interpreter' for Axes
    set(gca, 'TickDir', 'out', 'FontSize', 9, 'TickLabelInterpreter', 'none');
    
    % Overlay Values
    max_v = max(mat(:), [], 'omitnan');
    min_v = min(mat(:), [], 'omitnan');
    for i = 1:size(mat, 1)
        for j = 1:size(mat, 2)
            val = mat(i, j);
            if ~isnan(val)
                % Dynamic text color
                norm_v = (val - min_v) / (max_v - min_v + eps);
                if norm_v < 0.25 || norm_v > 0.8, col='w'; else, col='k'; end
                
                % Format text
                if val >= 10
                    txt = sprintf('%.1f', val); 
                else
                    txt = sprintf('%.2f', val); 
                end
                text(j, i, txt, 'HorizontalAlignment', 'center', ...
                     'Color', col, 'FontSize', 9, 'FontWeight', 'bold');
            end
        end
    end
end