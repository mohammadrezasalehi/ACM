%% Suggested Name: compare_real_vs_gaussian_null.m
% Compares Aggregated Network-Level Features: Real Brain vs. Gaussian Null.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';

% Paths to Network-Level Aggregated Features (All Subjects)
real_net_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features', 'EdgeInNetwork_Level');
gauss_net_dir = fullfile(project_root, 'results', 'Modulated_FC_Gaussian', 'Aggregated_Features', 'EdgeInNetwork_Level');

% Load Data
fprintf('Loading Real and Gaussian Data...\n');
RealData = load(fullfile(real_net_dir, 'AllSubjs_Net_Features.mat'));
GaussData = load(fullfile(gauss_net_dir, 'AllSubjs_Net_Features.mat'));

% We will average across all subjects and all sessions to get the ultimate Group Pattern
% RealData.feat_net_slope is an [100 x 4] cell array.
n_subs = size(RealData.feat_net_slope, 1);
n_sess = size(RealData.feat_net_slope, 2);

%% 2. Function to Extract Grand Mean Matrix [17 x 17]
get_grand_mean = @(cell_array) mean(cat(3, cell_array{:}), 3, 'omitnan');

% Extract Grand Means (Real)
R_var   = get_grand_mean(RealData.feat_net_var);
R_slope = get_grand_mean(RealData.feat_net_slope);
R_zcr   = get_grand_mean(RealData.feat_net_zcr);
R_poly  = get_grand_mean(RealData.feat_net_poly);

% Extract Grand Means (Gaussian)
G_var   = get_grand_mean(GaussData.feat_net_var);
G_slope = get_grand_mean(GaussData.feat_net_slope);
G_zcr   = get_grand_mean(GaussData.feat_net_zcr);
G_poly  = get_grand_mean(GaussData.feat_net_poly);

%% 3. Visualization: Scatter Plots (Real vs Gaussian)
figure('Name', 'Real vs Gaussian Null', 'Position', [50, 50, 1200, 800], 'Color', 'w');
sgtitle('Static FC Preserved (Gaussian) vs. True Non-Linear Dynamics (Real)', 'FontSize', 16, 'FontWeight', 'bold');

plot_comparison_scatter(1, G_var(:), R_var(:), 'Variance Sensitivity');
plot_comparison_scatter(2, G_slope(:), R_slope(:), 'Slope (Absolute Mean)');
plot_comparison_scatter(3, G_zcr(:), R_zcr(:), 'MCR/ZCR Volatility');
plot_comparison_scatter(4, G_poly(:), R_poly(:), 'Non-Linear Prevalence (%)');

%% --- HELPER FUNCTION ---
function plot_comparison_scatter(idx, gauss_vec, real_vec, tit)
    subplot(2, 2, idx);
    
    % Plot scatter points
    scatter(gauss_vec, real_vec, 40, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Calculate correlation
    r_val = corr(gauss_vec, real_vec, 'Rows', 'complete');
    
    % Fit a line to see if Real > Gaussian
    p = polyfit(gauss_vec, real_vec, 1);
    x_fit = linspace(min(gauss_vec), max(gauss_vec), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'b-', 'LineWidth', 2);
    
    % Plot Identity Line (y = x) in Red
    % If points are above the red line, the Real brain has STRONGER effects than Gaussian
    min_all = min([gauss_vec; real_vec]);
    max_all = max([gauss_vec; real_vec]);
    plot([min_all, max_all], [min_all, max_all], 'r--', 'LineWidth', 2);
    
    title(sprintf('%s\n(r = %.2f)', tit, r_val), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Gaussian Null (Preserves Static FC)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Real Brain Dynamics', 'FontSize', 10, 'FontWeight', 'bold');
    
    legend('Network Edges', 'Linear Fit', 'Identity Line (y=x)', 'Location', 'best');
    grid on;
    set(gca, 'TickDir', 'out');
    hold off;
end