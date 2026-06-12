% Generates a Bubble Chart of Node-Level Modulatory Features.
% X: Variance, Y: Non-Linearity, Size: ZCR/MCR, Color: Resting-State Network

clc; clear; close all;

%% 1. Configuration & Load Data
project_root = 'F:\PhD Code\My_PhD_Project';
base_result_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features');

% Load Node Features (Averaged across 4 sessions)
NodeData = load(fullfile(base_result_dir, 'Node_Level', 'AllSubjs_Node_Features.mat'));

% Calculate Grand Mean across all 100 subjects for each node (1x200)
var_mean  = mean(mean(cat(3, NodeData.feat_node_var{:}), 3), 1)';
poly_mean = mean(mean(cat(3, NodeData.feat_node_poly{:}), 3), 1)';
zcr_mean  = mean(mean(cat(3, NodeData.feat_node_zcr{:}), 3), 1)';

% Load Atlas Info (for coloring the 200 nodes)
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);

% The 7 Broad Yeo Networks we want to map to
unique_7_nets = {'Visual', 'Somatomotor', 'Dorsal Attention', 'Ventral Attention', 'Limbic', 'Control', 'Default Mode'};
n_broad = length(unique_7_nets);

net_indices = zeros(200, 1);

for i = 1:200
    net_name = roi_nets{i};
    
    % Map the 17 network names to the 7 broad categories safely
    if contains(net_name, 'Vis', 'IgnoreCase', true)
        net_indices(i) = 1; % Visual
    elseif contains(net_name, 'SomMot', 'IgnoreCase', true)
        net_indices(i) = 2; % Somatomotor
    elseif contains(net_name, 'DorsAttn', 'IgnoreCase', true)
        net_indices(i) = 3; % Dorsal Attention
    elseif contains(net_name, 'SalVenAttn', 'IgnoreCase', true)
        net_indices(i) = 4; % Ventral Attention
    elseif contains(net_name, 'Lim', 'IgnoreCase', true)
        net_indices(i) = 5; % Limbic
    elseif contains(net_name, 'Cont', 'IgnoreCase', true)
        net_indices(i) = 6; % Control
    elseif contains(net_name, 'Default', 'IgnoreCase', true)
        net_indices(i) = 7; % Default Mode
    else
        % Fallback for any unknown labels (shouldn't happen with Schaefer atlas)
        warning('Unrecognized network at node %d: %s', i, net_name);
        net_indices(i) = 6; % Default fallback to Control
    end
end

%% 2. Visualization Settings
% Standard Yeo 7-Network Colors
yeo_colors = [
    0.47, 0.18, 0.56; % Visual (Purple)
    0.27, 0.51, 0.71; % Somatomotor (Blue)
    0.00, 0.46, 0.11; % Dorsal Attention (Green)
    0.77, 0.23, 0.81; % Ventral Attention (Violet/Pink)
    0.86, 0.97, 0.64; % Limbic (Cream)
    0.90, 0.58, 0.13; % Control (Orange)
    0.80, 0.24, 0.27  % Default Mode (Red)
];

% Scale ZCR for Bubble Sizes (e.g., Min size = 30, Max size = 400)
min_size = 30; max_size = 400;
zcr_norm = (zcr_mean - min(zcr_mean)) / (max(zcr_mean) - min(zcr_mean));
bubble_sizes = min_size + zcr_norm * (max_size - min_size);

%% 3. Plotting the Bubble Chart
figure('Color', 'w', 'Position', [100, 100, 900, 700]);

hold on;
scatter_handles = zeros(n_broad, 1);

% Plot each network separately to build a clean legend
for n = 1:n_broad
    idx = (net_indices == n);
    
    % Use scatter with dynamic sizes and specified colors
    scatter_handles(n) = scatter(var_mean(idx), poly_mean(idx), bubble_sizes(idx), ...
        'MarkerFaceColor', yeo_colors(n, :), ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 0.75);
end

% Aesthetics
grid on;
box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1.2);

xlabel('Global Modulatory Strength (Node-Level Variance)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Non-Linear Footprint (Prevalence of Degree 2 & 3)', 'FontSize', 14, 'FontWeight', 'bold');
title('Modulatory Feature Space Across 200 Cortical Nodes', 'FontSize', 16, 'FontWeight', 'bold');

% Legend for Colors (Networks)
leg1 = legend(scatter_handles, unique_7_nets, 'Location', 'best', 'FontSize', 11);
title(leg1, 'Resting-State Networks');

% Optional: Add a text box explaining the bubble size
annotation('textbox', [0.15 0.8 0.3 0.1], 'String', 'Bubble Size = Volatility (MCR)', ...
    'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);

hold off;