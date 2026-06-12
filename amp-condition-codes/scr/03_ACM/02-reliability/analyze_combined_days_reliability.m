% Combines REST1 (Day 1) and REST2 (Day 2) by averaging extracted features.
% Calculates Day 1 vs. Day 2 Reliability (2x2 Heatmaps) for Node and Network levels.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
base_result_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features');

% Load Previously Extracted Features
fprintf('Loading previously extracted features...\n');
NodeData = load(fullfile(base_result_dir, 'Node_Level', 'AllSubjs_Node_Features.mat'));
NetData  = load(fullfile(base_result_dir, 'EdgeInNetwork_Level', 'AllSubjs_Net_Features.mat'));

[n_subs, n_sess] = size(NodeData.feat_node_var);

if n_sess ~= 4
    error('Expected 4 sessions (R1_LR, R1_RL, R2_LR, R2_RL). Check data.');
end

%% 2. Function to Combine Days and Calculate 2x2 Reliability
% Indices: 1=R1_LR, 2=R1_RL, 3=R2_LR, 4=R2_RL
% Day 1 = Mean of (1, 2)
% Day 2 = Mean of (3, 4)

combine_and_corr = @(feat_cell) compute_day_reliability(feat_cell, n_subs);

fprintf('Calculating Day 1 vs Day 2 Reliability...\n');
% --- Node Level ---
rel_node_var   = combine_and_corr(NodeData.feat_node_var);
rel_node_slope = combine_and_corr(NodeData.feat_node_slope);
rel_node_zcr   = combine_and_corr(NodeData.feat_node_zcr);
rel_node_poly  = combine_and_corr(NodeData.feat_node_poly);

% --- Network Level ---
rel_net_var    = combine_and_corr(NetData.feat_net_var);
rel_net_slope  = combine_and_corr(NetData.feat_net_slope);
rel_net_zcr    = combine_and_corr(NetData.feat_net_zcr);
rel_net_poly   = combine_and_corr(NetData.feat_net_poly);

%% 3. Visualization (2x2 Heatmaps)
labels = {'Day 1 (REST1)', 'Day 2 (REST2)'};

% --- Figure 1: Node-Level Combined Reliability ---
figure('Name', 'Node-Level Combined Reliability', 'Position', [50, 50, 1000, 800], 'Color', 'w');
sgtitle('Day-to-Day Reliability: Node-Level Modulator Power [1 x 200]', 'FontSize', 16, 'FontWeight', 'bold');
plot_rel_heatmap(1, rel_node_var, 'Variance (Standard Deviation)', labels);
plot_rel_heatmap(2, rel_node_slope, 'Slope (Absolute Mean)', labels);
plot_rel_heatmap(3, rel_node_zcr, 'MCR/ZCR (Mean)', labels);
plot_rel_heatmap(4, rel_node_poly, 'Poly Degree (Non-Linear %)', labels);

% --- Figure 2: Edge-In-Network Combined Reliability ---
figure('Name', 'Edge-in-Network Combined Reliability', 'Position', [100, 100, 1000, 800], 'Color', 'w');
sgtitle('Day-to-Day Reliability: Target Network Sensitivity [17 x 17]', 'FontSize', 16, 'FontWeight', 'bold');
plot_rel_heatmap(1, rel_net_var, 'Variance (Mean Sensitivity)', labels);
plot_rel_heatmap(2, rel_net_slope, 'Slope (Absolute Mean Sensitivity)', labels);
plot_rel_heatmap(3, rel_net_zcr, 'MCR/ZCR (Mean Sensitivity)', labels);
plot_rel_heatmap(4, rel_net_poly, 'Poly Degree (Non-Linear %)', labels);


%% --- HELPER FUNCTIONS ---

function avg_rel = compute_day_reliability(feat_cell, n_subs)
    % Calculates a 2x2 reliability matrix across subjects
    rel_3d = zeros(2, 2, n_subs);
    
    for i = 1:n_subs
        % Extract the 4 sessions for subject i
        s1 = feat_cell{i, 1}(:);
        s2 = feat_cell{i, 2}(:);
        s3 = feat_cell{i, 3}(:);
        s4 = feat_cell{i, 4}(:);
        
        % Combine (Average) Day 1 and Day 2
        day1 = (s1 + s2) / 2;
        day2 = (s3 + s4) / 2;
        
        data_mat = [day1, day2];
        
        % Calculate 2x2 correlation matrix
        rel_3d(:, :, i) = corr(data_mat, 'Rows', 'complete');
    end
    
    % Average across subjects
    avg_rel = mean(rel_3d, 3, 'omitnan');
end

function plot_rel_heatmap(idx, mat, tit, labels)
    subplot(2, 2, idx);
    imagesc(mat);
    colormap('jet'); colorbar;
    
    % Set limits
    min_val = min(mat(:));
    caxis([min(0.2, min_val * 0.9), 1]); 
    
    % Formatting
    xticks(1:2); yticks(1:2);
    xticklabels(labels); yticklabels(labels);
    set(gca, 'TickDir', 'out', 'FontSize', 11);
    title(tit, 'FontSize', 13, 'FontWeight', 'bold');
    
    % Overlay Numbers
    for r = 1:2
        for c = 1:2
            val = mat(r, c);
            if val < (min_val + 1)/2, t_col = 'w'; else, t_col = 'k'; end
            text(c, r, sprintf('%.3f', val), 'HorizontalAlignment', 'center', ...
                 'Color', t_col, 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
end