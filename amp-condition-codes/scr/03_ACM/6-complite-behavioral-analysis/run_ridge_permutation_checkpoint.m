%% run_ridge_permutation_checkpoint.m
% Generates Null distributions (1000 permutations) with Auto-Checkpointing.
% Safe to stop (Ctrl+C) and resume anytime to let the PC cool down.

clc; clear; close all;

%% 1. Configuration & Loading Data
project_root = 'F:\PhD Code\My_PhD_Project';
save_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features');

% --- Load Empirical Models ---
empirical_file = fullfile(save_dir, 'Ridge_Prediction_Results_5Factors.mat');
fprintf('Loading Empirical Models from:\n%s\n', empirical_file);
load(empirical_file, 'models_node', 'models_net', 'factor_names'); 

% --- Load Behavioral Targets ---
y_all_file = fullfile(save_dir, 'y_all.mat');
fprintf('Loading Y_all Vector from:\n%s\n', y_all_file);
load(y_all_file, 'Y_all');

n_subs = size(Y_all, 1);
n_models = length(models_node);
n_factors = length(factor_names);

%% 2. Permutation & Checkpoint Setup
n_perms = 1000;
k_folds = 10;
lambda = 10; % Must match the empirical script
save_interval = 20; % Save progress every 20 permutations (~12 mins)

null_save_path = fullfile(save_dir, 'Ridge_Permutation_Nulls_InProgress.mat');

% --- Checkpoint Logic ---
if isfile(null_save_path)
    % Resume from where it left off
    load(null_save_path, 'R_Null_Node', 'R_Null_Net', 'last_completed_p');
    start_p = last_completed_p + 1;
    fprintf('\nFound checkpoint! Resuming from permutation %d...\n', start_p);
else
    % Start fresh
    R_Null_Node = zeros(n_factors, n_models, n_perms);
    R_Null_Net  = zeros(n_factors, n_models, n_perms);
    start_p = 1;
    fprintf('\nNo checkpoint found. Starting fresh from permutation 1...\n');
end

if start_p > n_perms
    fprintf('All %d permutations are already completed! Done.\n', n_perms);
    return;
end

%% 3. Run Permutation Loop
fprintf('======================================================\n');
fprintf('Running Permutations (Saving every %d iterations)\n', save_interval);
fprintf('Tip: Press Ctrl+C anytime to pause. Run again to resume.\n');
fprintf('======================================================\n');

tic;
for p = start_p:n_perms
    
    % Use 'p' as the seed so that each permutation has a unique but reproducible shuffle
    % This ensures that if we stop and resume, we don't repeat the exact same shuffle
    rng(p, 'twister'); 
    
    % Shuffle Behavioral Data across subjects
    shuffle_idx = randperm(n_subs);
    Y_shuffled = Y_all(shuffle_idx, :);
    
    % Generate exactly one K-fold partition for this permutation
    cv_indices = crossvalind('Kfold', n_subs, k_folds);
    
    y_preds_node = zeros(n_subs, n_models, n_factors);
    y_preds_net  = zeros(n_subs, n_models, n_factors);
    
    % K-Fold Cross Validation
    for k = 1:k_folds
        test_idx  = (cv_indices == k);
        train_idx = (cv_indices ~= k);
        
        for f = 1:n_factors
            y_tr = Y_shuffled(train_idx, f);
            
            for m = 1:n_models
                % --- NODE LEVEL ---
                X_tr_node = models_node(m).X(train_idx, :);
                X_te_node = models_node(m).X(test_idx, :);
                
                % Normalize
                mu = mean(X_tr_node, 1); sig = std(X_tr_node, 0, 1); sig(sig==0)=1;
                X_tr_n = (X_tr_node - mu) ./ sig; X_te_n = (X_te_node - mu) ./ sig;
                
                % Train & Predict
                B_node = ridge(y_tr, X_tr_n, lambda, 0);
                y_preds_node(test_idx, m, f) = B_node(1) + X_te_n * B_node(2:end);
                
                % --- NETWORK LEVEL ---
                X_tr_net = models_net(m).X(train_idx, :);
                X_te_net = models_net(m).X(test_idx, :);
                
                % Normalize
                mu = mean(X_tr_net, 1); sig = std(X_tr_net, 0, 1); sig(sig==0)=1;
                X_tr_n = (X_tr_net - mu) ./ sig; X_te_n = (X_te_net - mu) ./ sig;
                
                % Train & Predict
                B_net = ridge(y_tr, X_tr_n, lambda, 0);
                y_preds_net(test_idx, m, f) = B_net(1) + X_te_n * B_net(2:end);
            end
        end
    end
    
    % Evaluate Predictions for this permutation
    for f = 1:n_factors
        y_true = Y_shuffled(:, f);
        for m = 1:n_models
            r_node = corr(y_preds_node(:, m, f), y_true, 'Type', 'Pearson');
            r_net  = corr(y_preds_net(:, m, f),  y_true, 'Type', 'Pearson');
            
            if isnan(r_node), r_node = 0; end
            if isnan(r_net),  r_net = 0;  end
            
            R_Null_Node(f, m, p) = r_node;
            R_Null_Net(f, m, p)  = r_net;
        end
    end
    
    % --- Auto-Save (Checkpointing) ---
    last_completed_p = p;
    if mod(p, save_interval) == 0 || p == n_perms
        save(null_save_path, 'R_Null_Node', 'R_Null_Net', 'last_completed_p', 'n_perms', '-v7.3');
        elapsed_time = toc;
        fprintf('Saved Progress: %d / %d permutations (%.1f%%) - Run time: %.1f min\n', ...
            p, n_perms, (p/n_perms)*100, elapsed_time/60);
    end
end

%% 4. Final Cleanup (Rename file to final version)
if last_completed_p == n_perms
    final_save_path = fullfile(save_dir, 'Ridge_Permutation_Nulls_1000.mat');
    movefile(null_save_path, final_save_path);
    fprintf('\nAll %d permutations completed successfully!\n', n_perms);
    fprintf('Final null distributions saved to:\n%s\n', final_save_path);
end