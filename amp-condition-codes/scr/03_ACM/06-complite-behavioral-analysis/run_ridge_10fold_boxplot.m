% Runs 10-Fold CV using Ridge Regression instead of CPM.
% Perfect for Node-Level (200) and Network-Level (289) features.

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


%% 4. Define Models (Same as before)
models_node = struct('name', {'1. Static', '2. Variance', '3. Slope', '4. ZCR', '5. PolyDegree', '6. All_Dyn'}, ...
    'X', {X_Node_Stat, X_Node_Var, X_Node_Slop, X_Node_ZCR, X_Node_Poly, [X_Node_Var, X_Node_Slop, X_Node_ZCR, X_Node_Poly]});

models_net = struct('name', {'1. Static', '2. Variance', '3. Slope', '4. ZCR', '5. PolyDegree', '6. All_Dyn'}, ...
    'X', {X_Net_Stat, X_Net_Var, X_Net_Slop, X_Net_ZCR, X_Net_Poly, [X_Net_Var, X_Net_Slop, X_Net_ZCR, X_Net_Poly]});

n_models = length(models_node);
factor_names = {'WellBeing', 'Cognition', 'Internalizing', 'ProcSpeed', 'SubstanceUse'};

%% 5. Run 10-Fold CV with Ridge Regression
k_folds = 10;
n_iterations = 20;

fprintf('\nRunning Ridge Regression for Node-Level Predictors...\n');
R_Node_Dist = run_repeated_ridge(models_node, Y_all, k_folds, n_iterations, factor_names);

fprintf('\nRunning Ridge Regression for Network-Level Predictors...\n');
R_Net_Dist  = run_repeated_ridge(models_net,  Y_all, k_folds, n_iterations, factor_names);

%% 6. Save Data for Statistical Testing and Plotting
fprintf('\nSaving Machine Learning Results...\n');
save_path = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features', 'Ridge_Prediction_Results_5Factors.mat');

% Pack everything needed for the next script
save(save_path, 'R_Node_Dist', 'R_Net_Dist', 'models_node', 'models_net', 'factor_names');

fprintf('Results saved successfully to:\n%s\n', save_path);

%% 7. Visualization
plot_cpm_boxplots(R_Node_Dist, models_node, factor_names, 'Node-Level Predictors (Ridge Regression)');
plot_cpm_boxplots(R_Net_Dist,  models_net,  factor_names, 'Edge-in-Network Predictors (Ridge Regression)');

%% ================= HELPER FUNCTIONS ================= %%

function R_Dist = run_repeated_ridge(models, Y_all, k_folds, n_iter, factor_names)
    n_subs = size(Y_all, 1);
    n_models = length(models);
    R_Dist = zeros(5, n_models, n_iter);
    
    for f = 1:5
        y_target = Y_all(:, f);
        
        for iter = 1:n_iter
            cv_indices = crossvalind('Kfold', n_subs, k_folds);
            y_preds = zeros(n_subs, n_models);
            
            for k = 1:k_folds
                test_idx  = (cv_indices == k);
                train_idx = (cv_indices ~= k);
                
                y_tr = y_target(train_idx);
                y_te = y_target(test_idx); % (??? ???? ??????? ????? ?? ????? ???? ??????? ??????)
                
                for m = 1:n_models
                    X = models(m).X;
                    X_tr = X(train_idx, :);
                    X_te = X(test_idx, :);
                    
                    % 1. Standardize Features (Crucial for Ridge!)
                    mu_X = mean(X_tr, 1);
                    sig_X = std(X_tr, 0, 1);
                    sig_X(sig_X == 0) = 1; % Prevent div by zero
                    
                    X_tr_norm = (X_tr - mu_X) ./ sig_X;
                    X_te_norm = (X_te - mu_X) ./ sig_X;
                    
                    % 2. Train Ridge Regression
                    % Using a fixed Lambda (e.g., 10 or 100) to prevent overfitting.
                    % In a full pipeline, Lambda can be tuned via nested CV.
                    lambda = 10; 
                    B = ridge(y_tr, X_tr_norm, lambda, 0); 
                    
                    % 3. Predict
                    % B(1) is the intercept, B(2:end) are the weights
                    y_preds(test_idx, m) = B(1) + X_te_norm * B(2:end);
                end
            end
            
            % Evaluate Iteration
            for m = 1:n_models
                r = corr(y_preds(:, m), y_target, 'Type', 'Pearson');
                if isnan(r), r = 0; end 
                R_Dist(f, m, iter) = r;
            end
        end
        fprintf('  -> %s done.\n', factor_names{f});
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