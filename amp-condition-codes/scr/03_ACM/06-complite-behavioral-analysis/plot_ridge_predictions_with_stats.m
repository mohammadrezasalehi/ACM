%% Suggested Name: plot_ridge_predictions_with_stats.m
% Loads the output of the Ridge Regression (10-Fold CV, 20 Iters).
% Calculates parametric P-values, applies FDR, and plots Boxplots with asterisks.

clc; clear; close all;

%% 1. Configuration & Load Data
project_root = 'F:\PhD Code\My_PhD_Project';
results_file = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features', 'Ridge_Prediction_Results_5Factors.mat');

fprintf('Loading Machine Learning Results...\n');
load(results_file); % Loads: R_Node_Dist, R_Net_Dist, models_node, models_net, factor_names

n_subs = 100; % Ensure this matches your actual number of subjects
n_factors = length(factor_names); % Should be 5
n_models = length(models_node);   % Should be 6

q_thresh = 0.05; % FDR threshold

%% 2. Process NODE-LEVEL Data
fprintf('Processing Node-Level Statistics...\n');
% Calculate Mean Accuracy across 20 iterations
Mean_R_Node = mean(R_Node_Dist, 3); % [5 x 6] matrix

% Calculate P-values and FDR
[P_Node, Sig_Node] = calculate_stats_and_fdr(Mean_R_Node, n_subs, q_thresh);

% Plot Node-Level Figure
plot_statistical_boxplots(R_Node_Dist, Mean_R_Node, Sig_Node, models_node, factor_names, ...
    'Behavioral Prediction (Node-Level Predictors [1 x 200])');


%% 3. Process NETWORK-LEVEL Data
fprintf('Processing Network-Level Statistics...\n');
% Calculate Mean Accuracy across 20 iterations
Mean_R_Net = mean(R_Net_Dist, 3); % [5 x 6] matrix

% Calculate P-values and FDR
[P_Net, Sig_Net] = calculate_stats_and_fdr(Mean_R_Net, n_subs, q_thresh);

% Plot Network-Level Figure
plot_statistical_boxplots(R_Net_Dist, Mean_R_Net, Sig_Net, models_net, factor_names, ...
    'Behavioral Prediction (Edge-in-Network Predictors [17 x 17])');

fprintf('Done!\n');


%% ================= HELPER FUNCTIONS ================= %%

function [P_mat, Sig_mat] = calculate_stats_and_fdr(Mean_R_mat, n_subs, q_thresh)
    [n_fac, n_mod] = size(Mean_R_mat);
    P_mat = ones(n_fac, n_mod); % Default to 1.0 (not significant)
    
    % 1. Calculate Parametric P-values
    df = n_subs - 2;
    for f = 1:n_fac
        for m = 1:n_mod
            r_val = Mean_R_mat(f, m);
            if r_val > 0
                % t-statistic for Pearson r
                t_stat = r_val * sqrt(df) / sqrt(1 - r_val^2 + eps);
                % 1-tailed p-value
                P_mat(f, m) = 1 - tcdf(t_stat, df);
            end
        end
    end
    
    % 2. Apply Benjamini-Hochberg FDR across ALL tests in the matrix (e.g., 5*6 = 30 tests)
    p_flat = P_mat(:);
    [p_sorted, sort_idx] = sort(p_flat);
    m_tests = length(p_sorted);
    fdr_thresholds = (1:m_tests)' / m_tests * q_thresh;
    
    sig_idx = find(p_sorted <= fdr_thresholds, 1, 'last');
    
    is_significant = false(size(p_flat));
    if ~isempty(sig_idx)
        is_significant(sort_idx(1:sig_idx)) = true;
    end
    
    % Reshape back to [5 x 6] matrix
    Sig_mat = reshape(is_significant, n_fac, n_mod);
end

function plot_statistical_boxplots(R_Dist, Mean_R, Sig_Matrix, models, factor_names, tit)
    figure('Color', 'w', 'Position', [50, 50, 1500, 500], 'Name', tit);
    set(gca, 'OuterPosition', [0,0,1,1])
    sgtitle(tit, 'FontSize', 16, 'FontWeight', 'bold');
    
    n_models = length(models);
    model_names = {models.name};
    n_factors = length(factor_names); % 5
    
    for f = 1:n_factors
        subplot(2, 5, f); % 2x3 grid accommodates 5 subplots
        
        % Data for this factor: Transpose [20 Iters x 6 Models]
        data_to_plot = squeeze(R_Dist(f, :, :))'; 
        
        % Draw Boxplot
        boxplot(data_to_plot, 'Labels', model_names, 'Colors', 'k', 'Symbol', 'o', 'Widths', 0.5);
        hold on;
        grid on; yline(0, 'k-', 'LineWidth', 1.5);
        
        % Overlay Mean Values and Asterisks
        for m = 1:n_models
            mean_val = Mean_R(f, m);
            
            % Draw Mean Value (Blue)
            text(m, max(data_to_plot(:, m)) + 0.04, sprintf('%.2f', mean_val), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
            
            % Draw Red Asterisk if significant
            if Sig_Matrix(f, m)
                text(m, max(data_to_plot(:, m)) + 0.08, '*', ...
                    'HorizontalAlignment', 'center', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'r');
            end
        end
        
        ax = findall(gcf, 'Type', 'axes');
        for k = 1:1
            pos = get(ax(k), 'Position');
            pos(2) = pos(2) - 0.05;   % move downward
%             pos(4) = pos(4) * 1;   % slightly reduce height
%             pos(1) = pos(1) + (f-2) * 0.05; %move back
            set(ax(k), 'Position', pos);
        end
        
                
        % Formatting
        xtickangle(30);
        if f==1
            ylabel('Pearson r (Accuracy)', 'FontWeight', 'bold', 'FontSize', 11);
        end
        title(factor_names{f}, 'FontSize', 14, 'FontWeight', 'bold');
        
        % Dynamic Y-Limits (Uniform but leaves space for stars)
        ylim([-0.25, max(0.4, max(data_to_plot(:))*1.3)]); 
        set(gca, 'TickDir', 'out', 'Box', 'off');
    end
end