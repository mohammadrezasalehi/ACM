% Extracts Node-Level and Edge-in-Network Level features for each subject.
% Calculates within-subject test-retest reliability across 4 sessions.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
sessions = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
clean_pipe = 'No_clean\'; % Adjust if necessary
base_result_dir = fullfile(project_root, 'results', 'Modulated_FC');

% Output Directories for Aggregated Features
out_dir_node = fullfile(base_result_dir, 'Aggregated_Features', 'Node_Level');
out_dir_net  = fullfile(base_result_dir, 'Aggregated_Features', 'EdgeInNetwork_Level');
if ~exist(out_dir_node, 'dir'), mkdir(out_dir_node); end
if ~exist(out_dir_net, 'dir'), mkdir(out_dir_net); end

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

% Network Index Map
net_idx_map = zeros(n_rois, 1);
for i = 1:n_rois
    net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i}));
end

% Get subject list from the first session
files = dir(fullfile(base_result_dir, sessions{1}, clean_pipe, '*_dyn_mod.mat'));
subject_names = {files.name};
n_subs = length(subject_names);
n_sess = length(sessions);

%% 2. Data Extraction & Aggregation Loop
% We will store the extracted vectors/matrices to avoid recalculating
fprintf('Extracting and Aggregating Features for %d subjects...\n', n_subs);

% Pre-allocate storage for reliability calculation [Subjects x Sessions]
% Using cell arrays to hold the matrices
feat_node_var   = cell(n_subs, n_sess);
feat_node_slope = cell(n_subs, n_sess);
feat_node_zcr   = cell(n_subs, n_sess);
feat_node_poly  = cell(n_subs, n_sess);

feat_net_var    = cell(n_subs, n_sess);
feat_net_slope  = cell(n_subs, n_sess);
feat_net_zcr    = cell(n_subs, n_sess);
feat_net_poly   = cell(n_subs, n_sess);

h = waitbar(0, 'Processing Subjects...');
for i = 1:n_subs
    subj_name = subject_names{i};
    
    for s = 1:n_sess
        % Load raw [19900 x 200] matrices
        file_path = fullfile(base_result_dir, sessions{s}, clean_pipe, subj_name);
        D = load(file_path);
        
        % Clean NaNs/255s temporarily for mathematical operations
        v_mat = D.variance; v_mat(isnan(v_mat)) = 0;
        t_mat = D.slope;    t_mat(isnan(t_mat)) = 0;
        z_mat = D.zcr;      z_mat(isnan(z_mat)) = 0;
        p_mat = D.best_poly_degree; 
        
        % ========================================================
        % PERSPECTIVE 1: NODE-LEVEL MODULATOR ROLE [1 x 200]
        % "How much does Node X control the rest of the brain?"
        % Average across edges (rows)
        % ========================================================
        feat_node_var{i, s}   = std(v_mat, 0, 1); % Using STD as discussed
        feat_node_slope{i, s} = mean(abs(t_mat), 1); % Absolute Mean
        feat_node_zcr{i, s}   = mean(z_mat, 1);
        
        % Poly Prevalence: % of edges with degree 2 or 3
        is_nonlin = (p_mat == 2 | p_mat == 3);
        valid_edges = sum(p_mat ~= 255, 1);
        feat_node_poly{i, s} = (sum(is_nonlin, 1) ./ max(valid_edges, 1)) * 100;
        
        % ========================================================
        % PERSPECTIVE 2: EDGE-IN-NETWORK TARGET SENSITIVITY [17 x 17]
        % "How much is the communication between Net A and Net B modulated?"
        % Average across modulators (columns), then group into 17x17
        % ========================================================
        v_mean_edge = mean(v_mat, 2);
        t_abs_mean_edge = mean(abs(t_mat), 2);
        z_mean_edge = mean(z_mat, 2);
        
        is_nonlin_edge = (p_mat == 2 | p_mat == 3);
        valid_mods = sum(p_mat ~= 255, 2);
        p_prev_edge = (sum(is_nonlin_edge, 2) ./ max(valid_mods, 1)) * 100;
        
        % Map [19900 x 1] to [17 x 17]
        feat_net_var{i, s}   = map_edges_to_network(v_mean_edge, u_map, v_map, net_idx_map, n_nets);
        feat_net_slope{i, s} = map_edges_to_network(t_abs_mean_edge, u_map, v_map, net_idx_map, n_nets);
        feat_net_zcr{i, s}   = map_edges_to_network(z_mean_edge, u_map, v_map, net_idx_map, n_nets);
        feat_net_poly{i, s}  = map_edges_to_network(p_prev_edge, u_map, v_map, net_idx_map, n_nets);
    end
    waitbar(i/n_subs, h, sprintf('Aggregating Data... %d%%', round((i/n_subs)*100)));
end
close(h);

% Optional: Save the extracted features for future use
save(fullfile(out_dir_node, 'AllSubjs_Node_Features.mat'), 'feat_node_var', 'feat_node_slope', 'feat_node_zcr', 'feat_node_poly');
save(fullfile(out_dir_net, 'AllSubjs_Net_Features.mat'), 'feat_net_var', 'feat_net_slope', 'feat_net_zcr', 'feat_net_poly');

%% 3. Calculate Subject-Level Reliability
fprintf('Calculating Test-Retest Reliability...\n');

% Function to compute average 4x4 reliability matrix across subjects
calc_rel = @(feat_cell) compute_avg_reliability(feat_cell, n_subs, n_sess);

rel_node_var = calc_rel(feat_node_var);
rel_node_slope = calc_rel(feat_node_slope);
rel_node_zcr = calc_rel(feat_node_zcr);
rel_node_poly = calc_rel(feat_node_poly);

rel_net_var = calc_rel(feat_net_var);
rel_net_slope = calc_rel(feat_net_slope);
rel_net_zcr = calc_rel(feat_net_zcr);
rel_net_poly = calc_rel(feat_net_poly);

%% 4. Visualization (Plotting Heatmaps)
sess_labels = {'R1\_LR', 'R1\_RL', 'R2\_LR', 'R2\_RL'};

% --- Figure 1: Node-Level Reliability ---
figure('Name', 'Node-Level Reliability', 'Position', [50, 50, 1000, 800], 'Color', 'w');
sgtitle('Subject Reliability: Node-Level Modulator Power [1 x 200]', 'FontSize', 16, 'FontWeight', 'bold');
plot_rel_heatmap(1, rel_node_var, 'Variance (Standard Deviation)', sess_labels);
plot_rel_heatmap(2, rel_node_slope, 'Slope (Absolute Mean)', sess_labels);
plot_rel_heatmap(3, rel_node_zcr, 'MCR/ZCR (Mean)', sess_labels);
plot_rel_heatmap(4, rel_node_poly, 'Poly Degree (Non-Linear %)', sess_labels);

% --- Figure 2: Edge-In-Network-Level Reliability ---
figure('Name', 'Edge-in-Network Reliability', 'Position', [100, 100, 1000, 800], 'Color', 'w');
sgtitle('Subject Reliability: Target Network Sensitivity [17 x 17]', 'FontSize', 16, 'FontWeight', 'bold');
plot_rel_heatmap(1, rel_net_var, 'Variance (Mean Sensitivity)', sess_labels);
plot_rel_heatmap(2, rel_net_slope, 'Slope (Absolute Mean Sensitivity)', sess_labels);
plot_rel_heatmap(3, rel_net_zcr, 'MCR/ZCR (Mean Sensitivity)', sess_labels);
plot_rel_heatmap(4, rel_net_poly, 'Poly Degree (Non-Linear %)', sess_labels);

fprintf('Analysis Complete!\n');

%% --- HELPER FUNCTIONS ---

function net_mat = map_edges_to_network(edge_vec, u, v, net_idx_map, n_nets)
    % Maps [19900x1] to [17x17] by averaging
    sum_scores = zeros(n_nets, n_nets);
    count_scores = zeros(n_nets, n_nets);
    
    for e = 1:length(u)
        ni = u(e); nj = v(e);
        net_i = net_idx_map(ni);
        net_j = net_idx_map(nj);
        val = edge_vec(e);
        
        sum_scores(net_i, net_j) = sum_scores(net_i, net_j) + val;
        count_scores(net_i, net_j) = count_scores(net_i, net_j) + 1;
        
        if net_i ~= net_j
            sum_scores(net_j, net_i) = sum_scores(net_j, net_i) + val;
            count_scores(net_j, net_i) = count_scores(net_j, net_i) + 1;
        end
    end
    net_mat = sum_scores ./ count_scores;
end

function avg_rel = compute_avg_reliability(feat_cell, n_subs, n_sess)
    rel_3d = zeros(n_sess, n_sess, n_subs);
    for i = 1:n_subs
        % Extract data for subject i across all 4 sessions
        % Flatten matrices into columns: [Num_Features x 4_Sessions]
        data_mat = zeros(numel(feat_cell{i, 1}), n_sess);
        for s = 1:n_sess
            data_mat(:, s) = feat_cell{i, s}(:);
        end
        % Calculate 4x4 correlation
        rel_3d(:, :, i) = corr(data_mat, 'Rows', 'complete');
    end
    % Average across subjects
    avg_rel = mean(rel_3d, 3, 'omitnan');
end

function plot_rel_heatmap(idx, mat, tit, labels)
    subplot(2, 2, idx);
    imagesc(mat);
    colormap('jet'); colorbar;
    
    % Dynamic caxis: Ensure max is 1, min is slightly below lowest value
    min_val = min(mat(:));
    caxis([min(0.2, min_val * 0.9), 1]); 
    
    % Formatting
    xticks(1:4); yticks(1:4);
    xticklabels(labels); yticklabels(labels);
    set(gca, 'TickDir', 'out', 'FontSize', 10);
    title(tit, 'FontSize', 12, 'FontWeight', 'bold');
    
    % Overlay Numbers
    for r = 1:4
        for c = 1:4
            val = mat(r, c);
            if val < (min_val + 1)/2, t_col = 'w'; else, t_col = 'k'; end
            text(c, r, sprintf('%.3f', val), 'HorizontalAlignment', 'center', ...
                 'Color', t_col, 'FontSize', 11, 'FontWeight', 'bold');
        end
    end
end