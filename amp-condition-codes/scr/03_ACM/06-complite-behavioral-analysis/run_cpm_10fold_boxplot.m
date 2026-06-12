% Runs 10-Fold Cross-Validated CPM for Static FC and Dynamic Features.
% Predicts 5 Latent Behavioral Factors.
% Generates Boxplots of prediction accuracy (Pearson r) over 20 iterations.

clc; clear; close all;

%% 1. Configuration & Paths
project_root = 'F:\PhD Code\My_PhD_Project';
sessions = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
clean_pipe = 'No_clean\'; % Adjust if your raw data is in another folder

base_result_dir = fullfile(project_root, 'results', 'Modulated_FC');
behav_file = fullfile(project_root, 'data', 'bihavioral', 'scores', 'scores_05.csv'); % Your behavior file

% Load Atlas Info (for Static FC Network Mapping)
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17; n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
n_edges = length(u_map);

net_idx_map = zeros(n_rois, 1);
for i = 1:n_rois
    net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i}));
end

%% 2. Load and Align Behavioral Data
fprintf('Aligning Behavioral Data...\n');
% Get Subject List from dynamic results
files = dir(fullfile(base_result_dir, sessions{1}, clean_pipe, '*_dyn_mod.mat'));
subject_names = {files.name};
n_subs = length(subject_names);
n_sess = length(sessions);

T_behav = readtable(behav_file);
factor_names = {'WellBeing', 'Cognition', 'Internalizing', 'ProcSpeed', 'SubstanceUse'};
Y_all = zeros(n_subs, 5);

for i = 1:n_subs
    % Extract Subject ID (assuming format like '100307_dyn_mod.mat' or 'sub-100307...')
    sid_str = regexp(subject_names{i}, '\d{6}', 'match', 'once'); 
    sid = str2double(sid_str);
    
    idx = find(T_behav.Subject == sid);
    if ~isempty(idx)
        Y_all(i, 1) = T_behav.factor1(idx);
        Y_all(i, 2) = T_behav.factor2(idx);
        Y_all(i, 3) = T_behav.factor3(idx);
        Y_all(i, 4) = T_behav.factor4(idx);
        Y_all(i, 5) = T_behav.factor5(idx);
    else
        warning('Subject %d not found in behavioral data! Filled with NaN.', sid);
        Y_all(i, :) = NaN;
    end
end

% Remove subjects with NaN behavior
valid_subs = ~isnan(Y_all(:,1));
Y_all = Y_all(valid_subs, :);
subject_names = subject_names(valid_subs);
n_subs = length(subject_names);

%% 3. Extract & Combine Features (4 Sessions Averaged)
fprintf('Aggregating 4-Session Features for %d subjects...\n', n_subs);

% Pre-allocate Feature Matrices [Subjects x Features]
% Node Level (Features = 200)
X_Node_Stat = zeros(n_subs, n_rois);
X_Node_Var  = zeros(n_subs, n_rois);
X_Node_Slop = zeros(n_subs, n_rois);
X_Node_ZCR  = zeros(n_subs, n_rois);
X_Node_Poly = zeros(n_subs, n_rois);

% Network Level (Features = 17x17 = 289)
n_net_feats = n_nets * n_nets;
X_Net_Stat = zeros(n_subs, n_net_feats);
X_Net_Var  = zeros(n_subs, n_net_feats);
X_Net_Slop = zeros(n_subs, n_net_feats);
X_Net_ZCR  = zeros(n_subs, n_net_feats);
X_Net_Poly = zeros(n_subs, n_net_feats);

% Load pre-extracted dynamic features
NodeData = load(fullfile(base_result_dir, 'Aggregated_Features', 'Node_Level', 'AllSubjs_Node_Features.mat'));
NetData  = load(fullfile(base_result_dir, 'Aggregated_Features', 'EdgeInNetwork_Level', 'AllSubjs_Net_Features.mat'));

h = waitbar(0, 'Calculating Static FC and Preparing Data...');
for i = 1:n_subs
    % --- Dynamic Features: Average across 4 sessions ---
    sub_orig_idx = find(valid_subs); % map back to 1-100 index
    real_i = sub_orig_idx(i);
    
    get_node = @(cell_mat) mean([cell_mat{real_i,1}(:), cell_mat{real_i,2}(:), cell_mat{real_i,3}(:), cell_mat{real_i,4}(:)], 2)';
    get_net  = @(cell_mat) mean([cell_mat{real_i,1}(:), cell_mat{real_i,2}(:), cell_mat{real_i,3}(:), cell_mat{real_i,4}(:)], 2)';
    
    X_Node_Var(i,:)  = get_node(NodeData.feat_node_var);
    X_Node_Slop(i,:) = get_node(NodeData.feat_node_slope);
    X_Node_ZCR(i,:)  = get_node(NodeData.feat_node_zcr);
    X_Node_Poly(i,:) = get_node(NodeData.feat_node_poly);
    
    X_Net_Var(i,:)  = get_net(NetData.feat_net_var);
    X_Net_Slop(i,:) = get_net(NetData.feat_net_slope);
    X_Net_ZCR(i,:)  = get_net(NetData.feat_net_zcr);
    X_Net_Poly(i,:) = get_net(NetData.feat_net_poly);
    
    % --- Static FC: Calculate on the fly for 4 sessions ---
    node_stat_tmp = zeros(4, n_rois);
    net_stat_tmp  = zeros(4, n_net_feats);
    
    for s = 1:n_sess
        % NOTE: Adjust this path if raw time-series are stored elsewhere
        raw_path = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', sessions{s}, subject_names{i});
        raw_path = strrep(raw_path, '_dyn_mod', ''); % Ensure correct raw filename
        tmp = load(raw_path);
        ts = tmp.data_on_atlas'; 
        
        static_fc = corr(ts);
        static_vec = static_fc(triu(true(n_rois), 1));
        
        % Node Level Static (Mean Absolute FC Strength)
        node_stat_tmp(s, :) = mean(abs(static_fc), 2)';
        % Network Level Static
        net_mat = map_edges_to_network(static_vec, u_map, v_map, net_idx_map, n_nets);
        net_stat_tmp(s, :) = net_mat(:)';
    end
    X_Node_Stat(i,:) = mean(node_stat_tmp, 1);
    X_Net_Stat(i,:)  = mean(net_stat_tmp, 1);
    
    waitbar(i/n_subs, h);
end
close(h);

%% 4. Define Models
models_node = struct('name', {'1. Static', '2. Variance', '3. Slope', '4. ZCR', '5. PolyDegree', '6. All_Dyn'}, ...
    'X', {X_Node_Stat, X_Node_Var, X_Node_Slop, X_Node_ZCR, X_Node_Poly, [X_Node_Var, X_Node_Slop, X_Node_ZCR, X_Node_Poly]});

models_net = struct('name', {'1. Static', '2. Variance', '3. Slope', '4. ZCR', '5. PolyDegree', '6. All_Dyn'}, ...
    'X', {X_Net_Stat, X_Net_Var, X_Net_Slop, X_Net_ZCR, X_Net_Poly, [X_Net_Var, X_Net_Slop, X_Net_ZCR, X_Net_Poly]});

n_models = length(models_node);
factor_names = {'WellBeing', 'Cognition', 'Internalizing', 'ProcSpeed', 'SubstanceUse'};

%% 5. Run 10-Fold CV CPM (Repeated)
k_folds = 10;
n_iterations = 20; % Repeat 10-fold CV 20 times for stability
p_thresh = 0.05;

fprintf('\nRunning 10-Fold CV for Node-Level Predictors...\n');
R_Node_Dist = run_repeated_cpm(models_node, Y_all, k_folds, n_iterations, p_thresh, factor_names);

fprintf('\nRunning 10-Fold CV for Network-Level Predictors...\n');
R_Net_Dist  = run_repeated_cpm(models_net,  Y_all, k_folds, n_iterations, p_thresh, factor_names);

%% 6. Visualization: Boxplots
plot_cpm_boxplots(R_Node_Dist, models_node, factor_names, 'Node-Level Predictors [1 x 200]');
plot_cpm_boxplots(R_Net_Dist,  models_net,  factor_names, 'Edge-in-Network Predictors [17 x 17]');

%% ================= HELPER FUNCTIONS ================= %%

function R_Dist = run_repeated_cpm(models, Y_all, k_folds, n_iter, p_thr, factor_names)
    n_subs = size(Y_all, 1);
    n_models = length(models);
    
    % Store correlation of each iteration: [5 Behaviors x 6 Models x 20 Iters]
    R_Dist = zeros(5, n_models, n_iter);
    
    for f = 1:5
        y_target = Y_all(:, f);
        
        for iter = 1:n_iter
            % Create random 10-fold partitions
            cv_indices = crossvalind('Kfold', n_subs, k_folds);
            
            y_preds = zeros(n_subs, n_models);
            
            for k = 1:k_folds
                test_idx  = (cv_indices == k);
                train_idx = (cv_indices ~= k);
                
                y_tr = y_target(train_idx);
                
                for m = 1:n_models
                    X = models(m).X;
                    y_preds(test_idx, m) = cpm_core_safe(X(train_idx, :), y_tr, X(test_idx, :), p_thr);
                end
            end
            
            % Evaluate Iteration (Correlation between predicted and actual)
            for m = 1:n_models
                r = corr(y_preds(:, m), y_target, 'Type', 'Pearson');
                % If CPM fails completely and predicts a flat line, r is NaN. Set to 0.
                if isnan(r), r = 0; end 
                R_Dist(f, m, iter) = r;
            end
        end
        fprintf('  -> %s done.\n', factor_names{f});
    end
end

function y_hat = cpm_core_safe(X_tr, y_tr, X_te, p_thr)
    % SAFE CPM: Prevents the -1.0 correlation bug
    [n_tr, ~] = size(X_tr);
    n_te = size(X_te, 1);
    
    r_vec = corr(X_tr, y_tr);
    r_vec(isnan(r_vec)) = 0;
    
    t_stat = r_vec .* sqrt(n_tr-2) ./ sqrt(1-r_vec.^2 + eps);
    p_vec = 2 * (1 - tcdf(abs(t_stat), n_tr-2));
    
    pos_mask = (r_vec > 0) & (p_vec < p_thr);
    neg_mask = (r_vec < 0) & (p_vec < p_thr);
    
    sum_tr = sum(X_tr(:, pos_mask), 2) - sum(X_tr(:, neg_mask), 2);
    sum_te = sum(X_te(:, pos_mask), 2) - sum(X_te(:, neg_mask), 2);
    
    % FIX: If no edges are selected, return exactly the mean of the training set for ALL test subjects.
    % This creates a flat prediction line -> Correlation will be exactly 0 (NaN handled later),
    % which correctly reflects that the model learned nothing, instead of giving a fake -1.0.
    if sum(pos_mask) == 0 && sum(neg_mask) == 0
        y_hat = repmat(mean(y_tr), n_te, 1);
    else
        p = polyfit(sum_tr, y_tr, 1);
        y_hat = polyval(p, sum_te);
    end
end

function plot_cpm_boxplots(R_Dist, models, factor_names, tit)
    figure('Color', 'w', 'Position', [50, 50, 1500, 800], 'Name', tit);
    sgtitle(sprintf('Behavioral Prediction (10-Fold CV, 20 Iters): %s', tit), 'FontSize', 16, 'FontWeight', 'bold');
    
    n_models = length(models);
    model_names = {models.name};
    
    for f = 1:5
        subplot(2, 3, f);
        
        % Extract data for this behavior: [6 Models x 20 Iters] -> Transpose to [20 x 6]
        data_to_plot = squeeze(R_Dist(f, :, :))'; 
        
        % Boxplot
        boxplot(data_to_plot, 'Labels', model_names, 'Colors', 'k', 'Symbol', 'o');
        hold on;
        
        % Add subtle grid and baseline
        grid on; 
        yline(0, 'r--', 'LineWidth', 1.5);
        
        % Find mean of each model to display as text above the box
        means = mean(data_to_plot, 1);
        for m = 1:n_models
            text(m, max(data_to_plot(:, m)) + 0.03, sprintf('%.2f', means(m)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
        end
        
        % Styling
        xtickangle(45);
        ylabel('Pearson r (Accuracy)', 'FontWeight', 'bold');
        title(factor_names{f}, 'FontSize', 14);
        
        % Uniform Y-limits for easy comparison
        ylim([-0.2, max(0.4, max(data_to_plot(:))*1.2)]);
        set(gca, 'TickDir', 'out');
    end
end

function net_mat = map_edges_to_network(edge_vec, u, v, net_idx_map, n_nets)
    sum_scores = zeros(n_nets, n_nets);
    count_scores = zeros(n_nets, n_nets);
    for e = 1:length(u)
        ni = u(e); nj = v(e);
        net_i = net_idx_map(ni); net_j = net_idx_map(nj);
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