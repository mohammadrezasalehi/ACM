%% run_ridge_combined_models.m
% Evaluates the incremental predictive value of Dynamic features 
% by combining them with Static FC (Static + Dynamic).

clc; clear; close all;

%% 1. Configuration & Loading Data
project_root = 'F:\PhD Code\My_PhD_Project';
save_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features');

% Load previously extracted features and behaviors
empirical_file = fullfile(save_dir, 'Ridge_Prediction_Results_5Factors.mat');
fprintf('Loading Base Features from:\n%s\n', empirical_file);
load(empirical_file, 'models_node', 'models_net', 'factor_names'); 

y_all_file = fullfile(save_dir, 'y_all.mat');
fprintf('Loading Behaviors from:\n%s\n', y_all_file);
load(y_all_file, 'Y_all');

n_subs = size(Y_all, 1);
n_factors = length(factor_names);

%% 2. Create Combined Models (Static + Dynamic)
fprintf('Building Combined (Concatenated) Feature Matrices...\n');

% Function to create combined struct
create_combined_struct = @(models) struct(...
    'name', {'1. Static', '2. Stat+Var', '3. Stat+Slp', '4. Stat+ZCR', '5. Stat+Poly', '6. Stat+All_Dyn'}, ...
    'X', { ...
        models(1).X, ...                                % 1. Static Only
        [models(1).X, models(2).X], ...                 % 2. Static + Variance
        [models(1).X, models(3).X], ...                 % 3. Static + Slope
        [models(1).X, models(4).X], ...                 % 4. Static + ZCR
        [models(1).X, models(5).X], ...                 % 5. Static + Poly
        [models(1).X, models(6).X]  ...                 % 6. Static + All Dynamic
    });

comb_models_node = create_combined_struct(models_node);
comb_models_net  = create_combined_struct(models_net);

%% 3. Run Ridge Regression for Combined Models
k_folds = 10;
n_iterations = 20;

fprintf('\nRunning Ridge CV for NODE-LEVEL Combined Models...\n');
R_Node_Comb = run_repeated_ridge(comb_models_node, Y_all, k_folds, n_iterations, factor_names);

fprintf('\nRunning Ridge CV for NETWORK-LEVEL Combined Models...\n');
R_Net_Comb  = run_repeated_ridge(comb_models_net,  Y_all, k_folds, n_iterations, factor_names);

%% 4. Plotting Results
fprintf('\nPlotting Boxplots...\n');
short_labels = {'St', 'St+Var', 'St+Slp', 'St+ZCR', 'St+Poly', 'St+All'};

plot_combined_boxplots(R_Node_Comb, short_labels, factor_names, ...
    'Incremental Value: Static + Dynamic (Node-Level [1 x 200])');

plot_combined_boxplots(R_Net_Comb, short_labels, factor_names, ...
    'Incremental Value: Static + Dynamic (Edge-in-Network [17 x 17])');

fprintf('Done!\n');

%% 5. Save Data for Statistical Testing and Plotting
fprintf('\nSaving Machine Learning Results...\n');
save_path = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features',...
    'Ridge_Prediction_Combined_Models.mat');

% Pack everything needed for the next script
save(save_path, 'R_Node_Comb', 'R_Net_Comb', 'comb_models_node', 'comb_models_net', 'Y_all', 'factor_names');

fprintf('Results saved successfully to:\n%s\n', save_path);


%% ================= HELPER FUNCTIONS ================= %%

function R_Dist = run_repeated_ridge(models, Y_all, k_folds, n_iter, factor_names)
    n_subs = size(Y_all, 1);
    n_models = length(models);
    R_Dist = zeros(5, n_models, n_iter);
    
    for f = 1:5
        y_target = Y_all(:, f);
        
        for iter = 1:n_iter
            % Fix seed based on iteration for reproducible folds
            rng(iter, 'twister'); 
            cv_indices = crossvalind('Kfold', n_subs, k_folds);
            y_preds = zeros(n_subs, n_models);
            
            for k = 1:k_folds
                test_idx  = (cv_indices == k);
                train_idx = (cv_indices ~= k);
                
                y_tr = y_target(train_idx);
                
                for m = 1:n_models
                    X = models(m).X;
                    X_tr = X(train_idx, :);
                    X_te = X(test_idx, :);
                    
                    % Normalize
                    mu_X = mean(X_tr, 1);
                    sig_X = std(X_tr, 0, 1);
                    sig_X(sig_X == 0) = 1; 
                    
                    X_tr_norm = (X_tr - mu_X) ./ sig_X;
                    X_te_norm = (X_te - mu_X) ./ sig_X;
                    
                    % Ridge
                    lambda = 10; 
                    B = ridge(y_tr, X_tr_norm, lambda, 0); 
                    
                    % Predict
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

function plot_combined_boxplots(R_Dist, x_labels, factor_names, tit)
    figure('Color', 'w', 'Position', [50, 50, 1500, 450], 'Name', tit);
    annotation('textbox', [0, 0.9, 1, 0.1], 'String', tit, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold');
    
    n_models = length(x_labels);
    n_factors = length(factor_names); 
    
    for f = 1:n_factors
        ax = subplot(1, 5, f); 
        
        data_to_plot = squeeze(R_Dist(f, :, :))'; 
        
        % Boxplot
        h = boxplot(ax, data_matrix_prep(data_to_plot), 'Labels', x_labels, 'Colors', [0.2 0.5 0.2], 'Symbol', 'o', 'Widths', 0.6);
        set(h, 'LineWidth', 1.2); 
        
        hold on;
        yline(0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        
        % Draw Mean Values
        means = mean(data_to_plot, 1);
        max_vals = max(data_to_plot, [], 1);
        offset = (max(data_to_plot(:)) - min(data_to_plot(:))) * 0.08;
        
        for m = 1:n_models
            text(m, max_vals(m) + offset, sprintf('%.2f', means(m)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k', 'FontWeight', 'bold');
        end
        
        % Aesthetics
        xtickangle(45);
        title(factor_names{f}, 'FontSize', 13, 'FontWeight', 'bold');
        set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
        
        if f == 1
            ylabel('Pearson r (Accuracy)', 'FontWeight', 'bold', 'FontSize', 11);
        end
        
        ylim([-0.25, max(0.4, max(max_vals) + offset*3)]);
    end
end

function out = data_matrix_prep(data)
    % Simple pass-through, just to ensure orientation is correct for boxplot
    out = data;
end