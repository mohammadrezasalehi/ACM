% Runs repeated 10-Fold CV Ridge Regression for predicting behavior.
% Computes parametric p-values for prediction accuracy (Pearson r).
% Corrects p-values using FDR and marks significant models with asterisks.

clc; clear; close all;

%% 1. Configuration & Data Load
% [Assuming X_Node_Stat, X_Net_Slop, Y_all, etc., are loaded exactly as before]
project_root = 'F:\PhD Code\My_PhD_Project';
factor_names = {'Mental Health', 'Cognition', 'ProcSpeed', 'Substance Use'};
n_factors = length(factor_names);

% Target models to evaluate
models_to_test = struct('name', {'1. Static', '2. Variance', '3. Slope', '4. ZCR', '5. PolyDegree', '6. All_Dyn'}, ...
    'X', {X_Net_Stat, X_Net_Var, X_Net_Slop, X_Net_ZCR, X_Net_Poly, [X_Net_Var, X_Net_Slop, X_Net_ZCR, X_Net_Poly]});
    
n_models = length(models_to_test);
n_subs = size(Y_all, 1);

%% 2. Run 10-Fold CV (Repeated)
k_folds = 10;
n_iterations = 20;

fprintf('1. Running Ridge Regression (10-Fold CV, %d Iters)...\n', n_iterations);
R_Dist_True = run_repeated_ridge(models_to_test, Y_all, k_folds, n_iterations);

% The "True" prediction accuracy is the MEAN accuracy across all iterations.
Mean_R_True = mean(R_Dist_True, 3); % Size: [n_factors x n_models]

%% 3. Parametric P-value Calculation
fprintf('2. Calculating Parametric P-values...\n');
P_parametric = zeros(n_factors, n_models);

for f = 1:n_factors
    for m = 1:n_models
        r_val = Mean_R_True(f, m);
        
        % If correlation is negative, it's not a successful prediction
        if r_val <= 0
            P_parametric(f, m) = 1.0;
        else
            % Calculate t-statistic for the Pearson correlation
            % Degrees of freedom = N - 2
            df = n_subs - 2;
            t_stat = r_val * sqrt(df) / sqrt(1 - r_val^2 + eps);
            
            % 1-tailed p-value (we only care if prediction is strictly POSITIVE)
            p_val = 1 - tcdf(t_stat, df);
            P_parametric(f, m) = p_val;
        end
    end
end

%% 4. FDR Correction (Benjamini-Hochberg)
fprintf('3. Applying FDR Correction...\n');
% Flatten P-values for global FDR correction across all tests (4 factors * 6 models = 24 tests)
p_flat = P_parametric(:);

% Manual FDR calculation
q_thresh = 0.05;
[p_sorted, sort_idx] = sort(p_flat);
m_tests = length(p_sorted);
fdr_thresholds = (1:m_tests)' / m_tests * q_thresh;
sig_idx = find(p_sorted <= fdr_thresholds, 1, 'last');

is_significant = false(size(p_flat));
if ~isempty(sig_idx)
    is_significant(sort_idx(1:sig_idx)) = true;
end

% Reshape back to [n_factors x n_models]
Significance_Matrix = reshape(is_significant, n_factors, n_models);

%% 5. Visualization (Boxplots with Significance Stars)
fprintf('4. Rendering Publication-Ready Figures...\n');
fig = figure('Color', 'w', 'Position', [50, 50, 1200, 800], 'Name', 'Ridge Predictions');
sgtitle('Behavioral Prediction Accuracy (Ridge Regression, Network-Level)', 'FontSize', 18, 'FontWeight', 'bold');

model_names = {models_to_test.name};

for f = 1:n_factors
    subplot(2, 2, f);
    
    % Transpose data to [Iters x Models] for boxplot
    data_to_plot = squeeze(R_Dist_True(f, :, :))'; 
    
    % Boxplot
    boxplot(data_to_plot, 'Labels', model_names, 'Colors', 'k', 'Symbol', 'o', 'Widths', 0.5);
    hold on;
    grid on; yline(0, 'k-', 'LineWidth', 1.5);
    
    % Plot Means and Stars
    for m = 1:n_models
        mean_val = Mean_R_True(f, m);
        
        % Draw Mean Value text (Blue)
        text(m, max(data_to_plot(:, m)) + 0.03, sprintf('%.2f', mean_val), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
        
        % Draw Red Asterisk if model survived FDR correction
        if Significance_Matrix(f, m)
            text(m, max(data_to_plot(:, m)) + 0.08, '*', ...
                'HorizontalAlignment', 'center', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'r');
        end
    end
    
    % Formatting
    xtickangle(45);
    ylabel('Pearson r (Predicted vs Actual)', 'FontWeight', 'bold', 'FontSize', 11);
    title(factor_names{f}, 'FontSize', 14, 'FontWeight', 'bold');
    
    % Dynamic Y-Limits
    ylim([-0.25, max(0.4, max(data_to_plot(:))*1.3)]); % Leave room for stars
    set(gca, 'TickDir', 'out', 'Box', 'off');
end

%% ================= HELPER FUNCTIONS ================= %%
function R_Dist = run_repeated_ridge(models, Y_all, k_folds, n_iter)
    n_subs = size(Y_all, 1);
    n_models = length(models);
    n_factors = size(Y_all, 2);
    R_Dist = zeros(n_factors, n_models, n_iter);
    
    for f = 1:n_factors
        y_target = Y_all(:, f);
        for iter = 1:n_iter
            cv_indices = crossvalind('Kfold', n_subs, k_folds);
            y_preds = zeros(n_subs, n_models);
            
            for k = 1:k_folds
                test_idx = (cv_indices == k); train_idx = ~test_idx;
                y_tr = y_target(train_idx);
                
                for m = 1:n_models
                    X = models(m).X;
                    X_tr = X(train_idx, :); X_te = X(test_idx, :);
                    
                    mu_X = mean(X_tr, 1); sig_X = std(X_tr, 0, 1); sig_X(sig_X == 0) = 1;
                    X_tr_norm = (X_tr - mu_X) ./ sig_X; X_te_norm = (X_te - mu_X) ./ sig_X;
                    
                    B = ridge(y_tr, X_tr_norm, 10, 0); 
                    y_preds(test_idx, m) = B(1) + X_te_norm * B(2:end);
                end
            end
            
            for m = 1:n_models
                r = corr(y_preds(:, m), y_target, 'Type', 'Pearson');
                if isnan(r), r = 0; end 
                R_Dist(f, m, iter) = r;
            end
        end
    end
end