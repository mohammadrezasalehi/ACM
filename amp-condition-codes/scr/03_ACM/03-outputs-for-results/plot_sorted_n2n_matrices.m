% Generates four 200x200 Node-by-Node Modulatory Matrices.
% Sorts nodes by the Kong 17-Network partition to reveal block topology.
% Adds network boundary lines and labels.

clc; clear; close all;

%% 1. Configuration & Load Data
project_root = 'F:\PhD Code\My_PhD_Project';
sessions = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
n_edges = length(u_map);

% Load Atlas Info for Sorting
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[roi_names, roi_nets] = load_atlas_info(atlas_info_path);

% Determine Network Membership and Unique Networks
unique_nets = unique(roi_nets, 'stable'); 
n_nets = length(unique_nets);

net_idx_map = zeros(n_rois, 1);
for i=1:n_rois, net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i})); end

% SORTING LOGIC: Sort nodes primarily by Network ID (1 to 17)
[sorted_net_idx, sort_order] = sort(net_idx_map);

% Find boundaries where network ID changes to draw lines later
net_boundaries = find(diff(sorted_net_idx) ~= 0);
boundary_ticks = [1; net_boundaries+1];

%% 2. Aggregate Data Across 4 Sessions
fprintf('Aggregating Group-Level Data into 200x200 matrices...\n');

Grand_Var = zeros(n_edges, n_rois);
Grand_Slop = zeros(n_edges, n_rois);
Grand_ZCR = zeros(n_edges, n_rois);
Grand_Poly = zeros(n_edges, n_rois); % Sum of Non-linear (Deg 2 & 3)

for s = 1:length(sessions)
    grp_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Group_Level', sessions{s});
    
    V = load(fullfile(grp_dir, 'Group_Variance_Consensus.mat'));  
    T = load(fullfile(grp_dir, 'Group_Slope_Tmap.mat'));          
    Z = load(fullfile(grp_dir, 'Group_ZCR_Median.mat'));          
    P = load(fullfile(grp_dir, 'Group_Poly_Prevalence.mat'));     
    
    Grand_Var  = Grand_Var + V.Consensus_Percent;
    Grand_Slop = Grand_Slop + T.T_map; % Retaining signs (+/-)
    Grand_ZCR  = Grand_ZCR + Z.Median;
    
    % Poly: Add prevalence of degrees 2 and 3
    is_nonlin = (P.Mode == 2 | P.Mode == 3);
    Grand_Poly = Grand_Poly + double(is_nonlin);
end

% Average over 4 sessions
Grand_Var  = Grand_Var / 4;
Grand_Slop = Grand_Slop / 4;
Grand_ZCR  = Grand_ZCR / 4;
Grand_Poly = (Grand_Poly / 4) * 100; % Convert back to %

%% 3. Convert Edge Vectors back to 200x200 Matrices
fprintf('Converting to Node-by-Node sorted format...\n');

Mat_Var = build_n2n_matrix(Grand_Var, u_map, v_map, n_rois);
Mat_Slop = build_n2n_matrix(Grand_Slop, u_map, v_map, n_rois);
Mat_ZCR = build_n2n_matrix(Grand_ZCR, u_map, v_map, n_rois);
Mat_Poly = build_n2n_matrix(Grand_Poly, u_map, v_map, n_rois);

% Apply the Sorting Order to Rows (Targets) and Columns (Modulators)
Sort_Var  = Mat_Var(sort_order, sort_order);
Sort_Slop = Mat_Slop(sort_order, sort_order);
Sort_ZCR  = Mat_ZCR(sort_order, sort_order);
Sort_Poly = Mat_Poly(sort_order, sort_order);

%% 4. Visualization
% Custom Colormaps
n_c = 256;
c_bwr = [linspace(0,1,n_c/2)', linspace(0,1,n_c/2)', linspace(1,1,n_c/2)'; 
         linspace(1,1,n_c/2)', linspace(1,0,n_c/2)', linspace(1,0,n_c/2)'];

fig = figure('Name', 'Sorted 200x200 Modulatory Matrices', 'Position', [50, 50, 1400, 1000], 'Color', 'w');
sgtitle('High-Resolution Modulatory Connectomes Sorted by 17 Functional Networks', 'FontSize', 18, 'FontWeight', 'bold');

% Subplot 1: Variance
plot_sorted_matrix(1, Sort_Var, 'hot', 'Modulatory Sensitivity (Variance)', boundary_ticks, unique_nets, false);

% Subplot 2: Slope (Bipolar: Red = Facil, Blue = Inhib)
plot_sorted_matrix(2, Sort_Slop, c_bwr, 'Directional Control (Mean Slope T-Value)', boundary_ticks, unique_nets, true);

% Subplot 3: ZCR
plot_sorted_matrix(3, Sort_ZCR, 'parula', 'Volatility (Mean-Crossing Rate)', boundary_ticks, unique_nets, false);

% Subplot 4: Polynomial Degree
plot_sorted_matrix(4, Sort_Poly, 'hot', 'Non-Linear Footprint (Prevalence %)', boundary_ticks, unique_nets, false);


%% ================= HELPER FUNCTIONS ================= %%
function n2n_mat = build_n2n_matrix(agg_vec_matrix, u, v, n_rois)
    % Averages across modulators (columns) to find how much each edge (row)
    % is modulated globally, then folds it back into a 200x200 target matrix.
    edge_means = mean(agg_vec_matrix, 2, 'omitnan');
    n2n_mat = zeros(n_rois, n_rois);
    for e = 1:length(u)
        ni = u(e); nj = v(e);
        n2n_mat(ni, nj) = edge_means(e);
        n2n_mat(nj, ni) = edge_means(e); % Symmetric
    end
    % Fill diagonal with NaN for aesthetic reasons
    n2n_mat(1:n_rois+1:end) = NaN; 
end

function plot_sorted_matrix(idx, mat, cmap, tit, b_ticks, net_labels, is_bipolar)
    subplot(2, 2, idx);
    imagesc(mat, 'AlphaData', ~isnan(mat));
    colormap(gca, cmap);
    
    if is_bipolar
        max_val = max(abs(mat(:)));
        if max_val == 0, max_val = 0.1; end
        caxis([-max_val, max_val]);
    end
    colorbar;
    
    axis square;
    title(tit, 'FontSize', 14, 'FontWeight', 'bold');
    
    % Draw Network Boundary Lines
    hold on;
    for i = 1:length(b_ticks)
        xline(b_ticks(i)-0.5, 'k-', 'LineWidth', 0.5);
        yline(b_ticks(i)-0.5, 'k-', 'LineWidth', 0.5);
    end
    
    % Set Ticks and Labels to Network Names
    % Calculate midpoints of networks for label placement
    label_pos = zeros(length(b_ticks), 1);
    for i = 1:length(b_ticks)-1
        label_pos(i) = (b_ticks(i) + b_ticks(i+1)) / 2;
    end
    label_pos(end) = (b_ticks(end) + 200) / 2;
    
    xticks(label_pos); yticks(label_pos);
    xticklabels(net_labels); yticklabels(net_labels);
    xtickangle(45);
    set(gca, 'TickDir', 'out', 'TickLabelInterpreter', 'none', 'FontSize', 7);
end