%% Suggested Name: run_residual_group_aggregation.m
% 1. Averages the subject-level Residual Matrices across all subjects.
% 2. Aggregates to Node and Network Levels.
% 3. Plots the Stacked R^2 Bar Chart and the Residual Brain Topographies.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

resid_dir = fullfile(project_root, 'results', 'Modulated_FC_Residuals', session, clean_pipe);
r2_file   = fullfile(project_root, 'results', 'Modulated_FC_Residuals', 'R2_Summary_Statistics.mat');

% Atlas Info
surf_path = fullfile(project_root, 'data', 'surf_info', 'surfinfo_kong2022_with_boundaries.mat'); 
load(surf_path); % sl, sr, gl, gr

n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
n_edges = length(u_map);

%% 2. Aggregate Residuals Across Subjects
files = dir(fullfile(resid_dir, '*_residualized.mat'));
n_subs = length(files);

fprintf('Aggregating Residuals for %d subjects...\n', n_subs);

Grand_Var  = zeros(n_edges, n_rois, 'single');
Grand_Slop = zeros(n_edges, n_rois, 'single');
Grand_ZCR  = zeros(n_edges, n_rois, 'single');
Grand_Poly = zeros(n_edges, n_rois, 'single');

for i = 1:n_subs
    R = load(fullfile(files(i).folder, files(i).name));
    
    % Replace NaNs with 0 for summation
    v=R.variance; v(isnan(v))=0; Grand_Var = Grand_Var + v;
    t=R.slope;    t(isnan(t))=0; Grand_Slop = Grand_Slop + t;
    z=R.zcr;      z(isnan(z))=0; Grand_ZCR = Grand_ZCR + z;
    p=R.best_poly_degree; p(isnan(p))=0; Grand_Poly = Grand_Poly + p;
end

Grand_Var  = Grand_Var / n_subs;
Grand_Slop = Grand_Slop / n_subs;
Grand_ZCR  = Grand_ZCR / n_subs;
Grand_Poly = Grand_Poly / n_subs;

%% 3. Node-Level Features of Residuals
fprintf('Calculating Node-Level Residual Topographies...\n');

Resid_Node_Var  = std(Grand_Var, 0, 1, 'omitnan')';
Resid_Node_ZCR  = mean(Grand_ZCR, 1, 'omitnan')';
Resid_Node_Poly = mean(Grand_Poly, 1, 'omitnan')';

% For Slope, residualization changes interpretation. 
% We look at the MEAN residual direction, retaining sign.
Resid_Node_Slop_Polarity = mean(Grand_Slop, 1, 'omitnan')'; 

%% 4. Visualization 1: Stacked R^2 Bar Chart
R2 = load(r2_file);
metrics = {'variance', 'slope', 'zcr', 'best_poly_degree'};
mean_FC = zeros(1, 4); mean_AR = zeros(1, 4);

for m = 1:4
    mean_FC(m) = mean(R2.(metrics{m}).R2_FC, 'omitnan');
    mean_AR(m) = mean(R2.(metrics{m}).R2_Full, 'omitnan') - mean_FC(m);
end

figure('Color', 'w', 'Position', [100, 100, 700, 500], 'Name', 'Variance Explained');
b = bar([mean_FC', mean_AR'], 'stacked', 'EdgeColor', 'k', 'LineWidth', 1.2);
b(1).FaceColor = [0.2 0.6 0.8]; b(2).FaceColor = [0.8 0.4 0.2];

xticklabels({'Variance', 'Slope', 'MCR/ZCR', 'Poly Degree'});
ylabel('Variance Explained (R^2)');
title('Variance Explained by Static FC vs Temporal Autocorrelation');
legend({'Model 1 (Static FC Triplets)', 'Model 2 Addition (\DeltaR^2 from AR1)'});
ylim([0, 1]); grid on; set(gca, 'TickDir', 'out', 'Box', 'off');

%% 5. Visualization 2: Residual Brain Maps
cmap_hot = hot(256);
rdbu = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];

figure('Color', 'w', 'Position', [150, 150, 1000, 800], 'Name', 'Residual Brain Maps');
sgtitle('Residual Topographies (Pure Dynamics after Regression)', 'FontSize', 16, 'FontWeight', 'bold');

plot_row_on_surf(1, 4, Resid_Node_Var, sr, sl, gr, gl, 'Residual Variance', cmap_hot);
plot_row_on_surf(2, 4, Resid_Node_Slop_Polarity, sr, sl, gr, gl, 'Residual Polarity', rdbu);
plot_row_on_surf(3, 4, Resid_Node_ZCR, sr, sl, gr, gl, 'Residual ZCR', parula(256));
plot_row_on_surf(4, 4, Resid_Node_Poly, sr, sl, gr, gl, 'Residual Poly Degree', cmap_hot);

%% 6. Statistical Comparison: Original vs Residual Spatial Maps
fprintf('\n========================================================\n');
fprintf('Comparing Original vs. Pure Dynamic Residual Topographies\n');
fprintf('========================================================\n');

% Load Original (Unadjusted) Group Data for Comparison
orig_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Group_Level', session);

Orig_V = load(fullfile(orig_dir, 'Group_Variance_Consensus.mat'));
Orig_T = load(fullfile(orig_dir, 'Group_Slope_Tmap.mat'));
Orig_Z = load(fullfile(orig_dir, 'Group_ZCR_Median.mat'));
Orig_P = load(fullfile(orig_dir, 'Group_Poly_Prevalence.mat'));

% Aggregate Original Data to Edge-Level (19900 x 1)
Orig_Edge_Var  = mean(Orig_V.Consensus_Percent, 2, 'omitnan');
Orig_Edge_Slop = mean(Orig_T.T_map, 2, 'omitnan');
Orig_Edge_ZCR  = mean(Orig_Z.Median, 2, 'omitnan');
Orig_Edge_Poly = mean(Orig_P.Prevalence_Percent, 2, 'omitnan');

% Aggregate Residual Data to Edge-Level (19900 x 1)
Resid_Edge_Var  = mean(Grand_Var, 2, 'omitnan');
Resid_Edge_Slop = mean(Grand_Slop, 2, 'omitnan');
Resid_Edge_ZCR  = mean(Grand_ZCR, 2, 'omitnan');
Resid_Edge_Poly = mean(Grand_Poly, 2, 'omitnan');

% Helper inline function to calculate and print correlations
calc_and_print_corr = @(orig, resid, name) ...
    print_corrs(orig, resid, name);

% Run correlations
calc_and_print_corr(Orig_Edge_Var, Resid_Edge_Var, 'Variance');
calc_and_print_corr(Orig_Edge_Slop, Resid_Edge_Slop, 'Slope (Directional Control)');
calc_and_print_corr(Orig_Edge_ZCR, Resid_Edge_ZCR, 'Volatility (MCR/ZCR)');
calc_and_print_corr(Orig_Edge_Poly, Resid_Edge_Poly, 'Poly Degree (Complexity)');

fprintf('========================================================\n');


%% 7. Visualization 3: Journal-Ready Statistical Heatmaps (Panel B & C)
fprintf('\nGenerating Statistical Heatmaps for Figure 3...\n');

% -------------------------------------------------------------------------
% PREPARE DATA FOR PANEL B (Variance Explained - R^2)
% Level of Aggregation: Raw Edge-Level Triplet Regression (Averaged across 100 subjects)
% -------------------------------------------------------------------------
% Metrics order: Variance, Slope, ZCR, Poly Degree
R2_FC_mean   = mean_FC';          % Column 1: FC Only
R2_Full_mean = mean_FC' + mean_AR'; % Column 2: FC + AR1
Delta_R2     = mean_AR';          % Column 3: Delta R^2

Data_Panel_B = [R2_FC_mean, R2_Full_mean, Delta_R2] * 100; % Convert to Percentage (%)

% -------------------------------------------------------------------------
% PREPARE DATA FOR PANEL C (Residual Preservation - Pearson r)
% -------------------------------------------------------------------------
% We calculate correlations between Original and Residual at 3 spatial levels.
% (Assuming you have Orig_Raw, Resid_Raw, Orig_Node, Resid_Node, etc. from earlier steps).
% For demonstration, I will calculate them here based on the aggregated variables we have:

% 1. Raw Level (19900 edges x 200 nodes averaged across subjects)
% (Needs the full 19900x200 matrices. If not in memory, replace with pre-calculated values)
% Here we use the Edge-Level (19900x1) as a proxy if full 2D is too heavy, 
% but ideally this is the raw 19900x200 correlation.
r_raw_var = corr(Orig_Edge_Var, Resid_Edge_Var, 'Rows', 'complete');
r_raw_slp = corr(Orig_Edge_Slop, Resid_Edge_Slop, 'Rows', 'complete');
r_raw_zcr = corr(Orig_Edge_ZCR, Resid_Edge_ZCR, 'Rows', 'complete');
r_raw_pol = corr(Orig_Edge_Poly, Resid_Edge_Poly, 'Rows', 'complete');

% 2. Node Level (200x1)
% (Assuming Orig_Node_Var and Resid_Node_Var exist in your workspace from Section 3)
r_nod_var = corr(Orig_Node_Var, Resid_Node_Var, 'Rows', 'complete');
r_nod_slp = corr(Orig_Node_Slop_Polarity, Resid_Node_Slop_Polarity, 'Rows', 'complete');
r_nod_zcr = corr(Orig_Node_ZCR, Resid_Node_ZCR, 'Rows', 'complete');
r_nod_pol = corr(Orig_Node_Poly, Resid_Node_Poly, 'Rows', 'complete');

% 3. Network Level (17x17)
% (Assuming Orig_Net_Var and Resid_Net_Var exist)
r_net_var = corr(Orig_Net_Var(:), Resid_Net_Var(:), 'Rows', 'complete');
r_net_slp = corr(Orig_Net_Slop(:), Resid_Net_Slop(:), 'Rows', 'complete');
r_net_zcr = corr(Orig_Net_ZCR(:), Resid_Net_ZCR(:), 'Rows', 'complete');
r_net_pol = corr(Orig_Net_Poly(:), Resid_Net_Poly(:), 'Rows', 'complete');

Data_Panel_C = [
    r_raw_var, r_nod_var, r_net_var;
    r_raw_slp, r_nod_slp, r_net_slp;
    r_raw_zcr, r_nod_zcr, r_net_zcr;
    r_raw_pol, r_nod_pol, r_net_pol
];

% -------------------------------------------------------------------------
% PLOTTING FUNCTION
% -------------------------------------------------------------------------
fig3_heatmaps = figure('Color', 'w', 'Position', [100, 100, 1100, 450], 'Name', 'Statistical Heatmaps');

metric_labels = {'Variance', 'Slope', 'ZCR', 'Non-Linear'};

% --- Plot Panel B: Variance Explained ---
subplot(1, 2, 1);
% Custom Colormap (White to Deep Orange/Red)
cmap_R2 = [ones(256,1), linspace(1,0.2,256)', linspace(1,0.2,256)'];
imagesc(Data_Panel_B); colormap(gca, cmap_R2); caxis([0, 100]);

% Formatting
title({'Panel B: Variance Explained (%)', '\fontsize{10}\fontweight{normal}(Raw Triplet-Level Regression)'}, 'FontSize', 13, 'FontWeight', 'bold');
xticks(1:3); xticklabels({'FC Only', 'FC + AR1', '\DeltaR^2 (AR1)'});
yticks(1:4); yticklabels(metric_labels);
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'off', 'XAxisLocation', 'top');
axis square;

% Overlay Text
for r = 1:4
    for c = 1:3
        val = Data_Panel_B(r, c);
        % Make text white if background is too dark, else black
        if val > 60, t_col = 'w'; else, t_col = 'k'; end
        text(c, r, sprintf('%.1f%%', val), 'HorizontalAlignment', 'center', ...
             'Color', t_col, 'FontSize', 12, 'FontWeight', 'bold');
    end
end
% Add Colorbar
cb1 = colorbar('Location', 'southoutside'); 
cb1.Label.String = 'Variance Explained (%)';

% --- Plot Panel C: Residual Preservation ---
subplot(1, 2, 2);
% Custom Colormap (White to Deep Blue for Correlations)
cmap_Corr = [linspace(1,0.1,256)', linspace(1,0.4,256)', ones(256,1)];
imagesc(Data_Panel_C); colormap(gca, cmap_Corr); caxis([0, 1]);

% Formatting
title({'Panel C: Residual Preservation', '\fontsize{10}\fontweight{normal}(Pearson r: Original vs Residualized)'}, 'FontSize', 13, 'FontWeight', 'bold');
xticks(1:3); xticklabels({'Raw', 'Node', 'Network'});
yticks(1:4); yticklabels(metric_labels);
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'off', 'XAxisLocation', 'top');
axis square;

% Overlay Text
for r = 1:4
    for c = 1:3
        val = Data_Panel_C(r, c);
        if val > 0.6, t_col = 'w'; else, t_col = 'k'; end
        text(c, r, sprintf('%.2f', val), 'HorizontalAlignment', 'center', ...
             'Color', t_col, 'FontSize', 12, 'FontWeight', 'bold');
    end
end
% Add Colorbar
cb2 = colorbar('Location', 'southoutside'); 
cb2.Label.String = 'Pearson Correlation (r)';

% Final Adjustments
% Adjust positions to bring them closer together
pos1 = get(subplot(1,2,1), 'Position');
pos2 = get(subplot(1,2,2), 'Position');
set(subplot(1,2,1), 'Position', [0.1, pos1(2), 0.35, pos1(4)]);
set(subplot(1,2,2), 'Position', [0.55, pos2(2), 0.35, pos2(4)]);

%% --- ADD THIS TO YOUR EXISTING HELPER FUNCTIONS ---
function print_corrs(orig_vec, resid_vec, metric_name)
    % Filter out NaNs and Infs
    valid_idx = isfinite(orig_vec) & isfinite(resid_vec);
    
    if sum(valid_idx) < 100
        fprintf('%-30s : Not enough valid data points.\n', metric_name);
        return;
    end
    
    % Calculate Correlations
    [r_pearson, p_pearson]   = corr(orig_vec(valid_idx), resid_vec(valid_idx), 'Type', 'Pearson');
    [r_spearman, p_spearman] = corr(orig_vec(valid_idx), resid_vec(valid_idx), 'Type', 'Spearman');
    
    % Print Results
    fprintf('%-30s:\n', metric_name);
    fprintf('   Pearson r  = %.4f (p = %.4e)\n', r_pearson, p_pearson);
    fprintf('   Spearman ? = %.4f (p = %.4e)\n\n', r_spearman, p_spearman);
end

%% --- HELPER FUNCTIONS ---
function plot_row_on_surf(row_idx, total_rows, plot_vals, sr, sl, gr, gl, tit, cmap)
    cr = zeros(size(gr.cdata)); cl = zeros(size(gl.cdata));
    cr(gr.cdata ~= 0) = plot_vals(gr.cdata(gr.cdata ~= 0));
    cl(gl.cdata ~= 0) = plot_vals(gl.cdata(gl.cdata ~= 0));
    
    clims = [min(plot_vals), max(plot_vals)];
    if diff(clims) == 0, clims = [0 1]; end
    if contains(tit, 'Polarity'), clims = [-max(abs(plot_vals)), max(abs(plot_vals))]; end 
    
    is_gifti = isfield(sr, 'data');
    idx_base = (row_idx - 1) * 4;
    
    ax1 = subplot(total_rows, 4, idx_base + 1); plot_surf(sl, cl, is_gifti, -90, clims); colormap(ax1, cmap);
    text(ax1, -0.4, 0.5, tit, 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    ax2 = subplot(total_rows, 4, idx_base + 2); plot_surf(sl, cl, is_gifti, 90, clims); colormap(ax2, cmap);
    ax3 = subplot(total_rows, 4, idx_base + 3); plot_surf(sr, cr, is_gifti, 90, clims); colormap(ax3, cmap);
    ax4 = subplot(total_rows, 4, idx_base + 4); plot_surf(sr, cr, is_gifti, -90, clims); colormap(ax4, cmap);
end

function plot_surf(s, c, is_gifti, az, clims)
    if is_gifti, th = trisurf(s.data{2}.data+1, s.data{1}.data(:,1), s.data{1}.data(:,2), s.data{1}.data(:,3), c);
    else, th = trisurf(s.faces, s.vertices(:,1), s.vertices(:,2), s.vertices(:,3), c); end
    set(th, 'edgecolor', 'none'); axis image off; set(gca, 'clim', clims);
    view(az, 0); material dull; camlight headlight; lighting gouraud;
end