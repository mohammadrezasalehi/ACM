%% plot_ridge_predictions_publication_ready.m
% Loads Empirical and Permutation results.
% Plots a highly readable, publication-ready figure.
% Stat Tests: Permutation vs Chance & Paired Permutation vs Static.

clc; clear; close all;

%% 1. Configuration & Load Data
project_root = 'F:\PhD Code\My_PhD_Project';
save_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features');

% Load Empirical Results (20 iterations)
empirical_file = fullfile(save_dir, 'Ridge_Prediction_Results_5Factors.mat');
fprintf('Loading Empirical Results...\n');
load(empirical_file, 'R_Node_Dist', 'R_Net_Dist', 'models_node', 'factor_names');

% Load Permutation Nulls (1000 permutations)
null_file = fullfile(save_dir, 'Ridge_Permutation_Nulls_1000.mat');
fprintf('Loading Permutation Nulls...\n');
load(null_file, 'R_Null_Node', 'R_Null_Net', 'n_perms');

n_factors = length(factor_names); % 5
n_models = length(models_node);   % 6

% Define Short Labels for X-axis to fix the "crowded" issue
short_labels = {'St', 'Var', 'Slp', 'ZCR', 'Poly', 'All'};

%% 2. Calculate Statistics
fprintf('Calculating Statistics (Paired Permutation Tests)...\n');

% --- Node Level Stats ---
[P_Chance_Node, Sig_Chance_Node, Sig_Stat_Node] = ...
    compute_dual_stats_perm(R_Node_Dist, R_Null_Node, n_perms, n_factors, n_models);

% --- Network Level Stats ---
[P_Chance_Net, Sig_Chance_Net, Sig_Stat_Net] = ...
    compute_dual_stats_perm(R_Net_Dist, R_Null_Net, n_perms, n_factors, n_models);

%% 3. Plotting the Figure
fprintf('Plotting Figure...\n');
% Create a large, high-resolution figure window
fig = figure('Color', 'w', 'Position', [50, 50, 1600, 900], 'Name', 'Behavioral Prediction');

% --- Row 1: Node-Level Predictors ---
annotation('textbox', [0, 0.92, 1, 0.05], 'String', 'A. Behavioral Prediction (Node-Level Predictors [1 x 200])', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');

for f = 1:n_factors
    ax = subplot(2, 5, f);
    plot_single_boxplot(ax, squeeze(R_Node_Dist(f,:,:))', short_labels, factor_names{f}, ...
        Sig_Chance_Node(f,:), Sig_Stat_Node(f,:), f==1);
end

% --- Row 2: Edge-in-Network Predictors ---
annotation('textbox', [0, 0.45, 1, 0.05], 'String', 'B. Behavioral Prediction (Edge-in-Network Predictors [17 x 17])', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');

for f = 1:n_factors
    ax = subplot(2, 5, f + 5);
    plot_single_boxplot(ax, squeeze(R_Net_Dist(f,:,:))', short_labels, factor_names{f}, ...
        Sig_Chance_Net(f,:), Sig_Stat_Net(f,:), f==1);
end

% --- Add Unified Legend at the Bottom ---
add_unified_legend(fig);

fprintf('Done! Save the figure as PDF or high-res TIFF.\n');

%% ================= HELPER FUNCTIONS ================= %%

function [P_Chance, Sig_Chance, Sig_Static] = compute_dual_stats_perm(Emp_Dist, Null_Dist, n_perms, n_fac, n_mod)
    % Emp_Dist: [Factors x Models x 20 iterations]
    % Null_Dist: [Factors x Models x 1000 permutations]
    
    Mean_Emp = mean(Emp_Dist, 3);
    P_Chance = ones(n_fac, n_mod);
    Sig_Chance = false(n_fac, n_mod);
    Sig_Static = false(n_fac, n_mod);
    
    for f = 1:n_fac
        for m = 1:n_mod
            % 1. Vs Chance (Standard Permutation)
            emp_val = Mean_Emp(f, m);
            null_vals = squeeze(Null_Dist(f, m, :));
            p_chance = (1 + sum(null_vals >= emp_val)) / (1 + n_perms);
            P_Chance(f, m) = p_chance;
            Sig_Chance(f, m) = (p_chance < 0.05);
            
            % 2. Vs Static (Paired Permutation Test)
            % Implemented exactly based on your logic:
            if m > 1
                % Empirical Difference (Dynamic - Static)
                emp_diff = Mean_Emp(f, m) - Mean_Emp(f, 1);
                
                % Null Distribution of Differences (Since permutations are paired by seed/iteration)
                null_diff = squeeze(Null_Dist(f, m, :)) - squeeze(Null_Dist(f, 1, :));
                
                % Non-parametric P-value for the difference
                p_stat = (1 + sum(null_diff >= emp_diff)) / (1 + n_perms);
                Sig_Static(f, m) = (p_stat < 0.05);
            end
        end
    end
end

function plot_single_boxplot(ax, data_matrix, x_labels, title_str, sig_chance, sig_static, is_first_col)
    % Draw the boxplot
    h = boxplot(ax, data_matrix, 'Labels', x_labels, 'Colors', [0.3 0.3 0.3], 'Symbol', 'o', 'Widths', 0.6);
    set(h, 'LineWidth', 1.2); 
    
    hold on;
    yline(0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    
    n_mod = size(data_matrix, 2);
    means = mean(data_matrix, 1);
    
    max_vals = max(data_matrix, [], 1);
    y_range = max(data_matrix(:)) - min(data_matrix(:));
    offset = y_range * 0.08; 
    
    for m = 1:n_mod
        % Print Mean Value
        text(m, max_vals(m) + offset*0.5, sprintf('%.2f', means(m)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0 0.4 0.7], 'FontWeight', 'bold');
        
        % Add Significance vs Chance (*)
        if sig_chance(m)
            text(m, max_vals(m) + offset*1.3, '*', ...
                'HorizontalAlignment', 'center', 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
        end
        
        % Add Significance vs Static (+)
        if sig_static(m)
            text(m, max_vals(m) + offset*2.1, '+', ...
                'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', [0.8 0 0], 'FontWeight', 'bold');
        end
    end
    
    title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10, 'LineWidth', 1);
    
    if is_first_col
        ylabel('Pearson r (Accuracy)', 'FontWeight', 'bold', 'FontSize', 12);
    end
    
    upper_lim = max(0.4, max(max_vals) + offset*3.5);
    ylim([-0.25, upper_lim]);
end

function add_unified_legend(fig)
    lax = axes('Position', [0 0 1 0.1], 'Visible', 'off');
    hold(lax, 'on');
    
    p1 = plot(lax, [NaN NaN], 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
    p2 = plot(lax, NaN, NaN, '*k', 'MarkerSize', 10);
    p3 = plot(lax, NaN, NaN, '+r', 'MarkerSize', 8, 'LineWidth', 2);
    
    leg = legend([p1, p2, p3], ...
        {'Predictors: St=Static, Var=Variance, Slp=Slope, ZCR=Volatility, Poly=Non-Linear, All=Combined', ...
         '* Significant prediction (Permutation p < 0.05 vs. Chance)', ...
         '+ Significantly outperforms Static FC (Paired Permutation p < 0.05)'}, ...
         'Orientation', 'horizontal', 'FontSize', 11, 'Box', 'off');
     
    leg.ItemTokenSize = [15, 18];
end