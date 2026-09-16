%% Suggested Name: run_spatial_regression_and_residual_maps_v2.m
% 1. Uses Sufficient Statistics (X'X, X'Y) to fit exact OLS on ~4M points.
% 2. Fits Nested Models: M_FC and M_Full (FC + AR1).
% 3. Calculates exact R^2 for FC and Delta R^2 for AR1.
% 4. Generates Residual Topographies from the Full Model.

clc; clear; close all;

%% 1. Configuration & Load Data
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

pred_dir  = fullfile(project_root, 'results', 'Modulated_FC', 'Control_Predictors', session, clean_pipe);
dyn_dir   = fullfile(project_root, 'results', 'Modulated_FC', session, clean_pipe);
surf_path = fullfile(project_root, 'data', 'surf_info', 'surfinfo_kong2022_with_boundaries.mat'); 
load(surf_path); % sl, sr, gl, gr

fprintf('Loading Predictor Matrices...\n');
PredData = load(fullfile(pred_dir, 'AllSubjs_Static_Predictors_Triplet.mat'));
n_subs = length(PredData.FC_XY);
[n_edges, n_rois] = size(PredData.FC_XY{1});
n_elements = n_edges * n_rois; 

dyn_files = dir(fullfile(dyn_dir, '*_dyn_mod.mat'));

%% 2. Process Subjects (Nested Models & Exact OLS)
fprintf('Running Exact Spatial Regression via Sufficient Statistics...\n');

% Storage: [Subjects x 4 Metrics]
R2_FC_Values   = zeros(n_subs, 4);
R2_Full_Values = zeros(n_subs, 4);

% Accumulators for Grand Mean Residual Maps
Grand_Resid_Var  = zeros(n_edges, n_rois, 'single');
Grand_Resid_Slop = zeros(n_edges, n_rois, 'single');
Grand_Resid_ZCR  = zeros(n_edges, n_rois, 'single');
Grand_Resid_Poly = zeros(n_edges, n_rois, 'single');

h = waitbar(0, 'Regressing out Static features...');
for i = 1:n_subs
    % 1. Load Dynamic Data (Y)
    D = load(fullfile(dyn_files(i).folder, dyn_files(i).name));
    
    Y_mats = {D.variance(:), D.slope(:), D.zcr(:), single(D.best_poly_degree(:))};
    valid_mask = ~isnan(Y_mats{1}) & ~isnan(Y_mats{2}) & ~isnan(Y_mats{3}) & (Y_mats{4} ~= 255);
    n_valid = sum(valid_mask);
    
    % 2. Build Design Matrices (X)
    % Z-scoring predictors is crucial to avoid ill-conditioned X'X matrix
    X_FC = [ones(n_valid, 1), ...
            zscore(PredData.FC_XY{i}(valid_mask)), ...
            zscore(PredData.FC_XZ{i}(valid_mask)), ...
            zscore(PredData.FC_YZ{i}(valid_mask))]; % Model 1: Intercept + 3 FCs
            
    X_Full = [X_FC, ...
              zscore(PredData.AR_X{i}(valid_mask)), ...
              zscore(PredData.AR_Y{i}(valid_mask)), ...
              zscore(PredData.AR_Z{i}(valid_mask))]; % Model 2: Model 1 + 3 AR1s
    
    % 3. Fit Models for each Metric using Sufficient Statistics
    resid_mats = cell(1, 4);
    
    for m = 1:4
        Y_tgt = Y_mats{m}(valid_mask);
        
        % Calculate Total Sum of Squares (SST) for R^2
        SST = sum((Y_tgt - mean(Y_tgt)).^2);
        
        % --- MODEL 1: FC Only ---
        Beta_FC = (X_FC' * X_FC) \ (X_FC' * Y_tgt);
        Y_hat_FC = X_FC * Beta_FC;
        SS_res_FC = sum((Y_tgt - Y_hat_FC).^2);
        R2_FC_Values(i, m) = 1 - (SS_res_FC / SST);
        
        % --- MODEL 2: Full (FC + AR1) ---
        Beta_Full = (X_Full' * X_Full) \ (X_Full' * Y_tgt);
        Y_hat_Full = X_Full * Beta_Full;
        Resid_Full = Y_tgt - Y_hat_Full;
        SS_res_Full = sum(Resid_Full.^2);
        R2_Full_Values(i, m) = 1 - (SS_res_Full / SST);
        
        % Store Residuals from FULL model for brain maps
        Full_Resid = NaN(n_elements, 1, 'single');
        Full_Resid(valid_mask) = Resid_Full;
        resid_mats{m} = reshape(Full_Resid, n_edges, n_rois);
    end
    
    % Add to Grand Mean Accumulators
    Grand_Resid_Var  = Grand_Resid_Var  + resid_mats{1};
    Grand_Resid_Slop = Grand_Resid_Slop + resid_mats{2};
    Grand_Resid_ZCR  = Grand_Resid_ZCR  + resid_mats{3};
    Grand_Resid_Poly = Grand_Resid_Poly + resid_mats{4};
    
    waitbar(i/n_subs, h);
end
close(h);

Grand_Resid_Var  = Grand_Resid_Var / n_subs;
Grand_Resid_Slop = Grand_Resid_Slop / n_subs;
Grand_Resid_ZCR  = Grand_Resid_ZCR / n_subs;
Grand_Resid_Poly = Grand_Resid_Poly / n_subs;

%% 3. Node-Level Aggregation of Residuals
fprintf('Aggregating Residuals to Node-Level...\n');
Resid_Node_Var  = std(Grand_Resid_Var, 0, 1, 'omitnan')';
Resid_Node_Facil= sum(Grand_Resid_Slop > 2.5, 1)'; 
Resid_Node_Inhib= sum(Grand_Resid_Slop < -2.5, 1)';
Resid_Node_Polr = (Resid_Node_Facil - Resid_Node_Inhib) ./ (Resid_Node_Facil + Resid_Node_Inhib + eps);
Resid_Node_ZCR  = mean(Grand_Resid_ZCR, 1, 'omitnan')';
Resid_Node_Poly = mean(Grand_Resid_Poly, 1, 'omitnan')'; 

%% 4. Visualization 1: Nested R^2 Bar Chart (Stacked)
fig1 = figure('Color', 'w', 'Position', [100, 100, 800, 600], 'Name', 'Nested R^2');

mean_R2_FC = mean(R2_FC_Values, 1);
mean_R2_AR = mean(R2_Full_Values, 1) - mean_R2_FC; % Delta R^2

% Stacked Bar Chart to show Total R^2 and its components
b = bar([mean_R2_FC', mean_R2_AR'], 'stacked', 'EdgeColor', 'k', 'LineWidth', 1.2);
b(1).FaceColor = [0.2 0.6 0.8]; % Blue for FC
b(2).FaceColor = [0.8 0.4 0.2]; % Orange for AR1

metric_names = {'Variance', 'Slope', 'MCR/ZCR', 'Poly Degree'};
xticks(1:4); xticklabels(metric_names);
ylabel('Variance Explained (R^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Variance Explained by Static FC vs Temporal Autocorrelation', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model 1 (Static FC Triplets)', 'Model 2 Addition (\DeltaR^2 from AR1)'}, 'Location', 'NorthWest');
ylim([0, 1]); grid on; set(gca, 'TickDir', 'out', 'Box', 'off');

% Add text values for Total R^2
total_R2 = mean_R2_FC + mean_R2_AR;
for m = 1:4
    text(m, total_R2(m) + 0.03, sprintf('%.1f%%', total_R2(m)*100), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

%% 5. Visualization 2: Residual Brain Topographies
cmap_hot = hot(256);
rdbu = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];

fig2 = figure('Color', 'w', 'Position', [150, 150, 1000, 800], 'Name', 'Residual Brain Maps');
sgtitle('Residual Topographies (Pure Dynamics after full regression)', 'FontSize', 16, 'FontWeight', 'bold');

plot_row_on_surf(1, 4, Resid_Node_Var, sr, sl, gr, gl, 'Residual Variance', cmap_hot);
plot_row_on_surf(2, 4, Resid_Node_Polr, sr, sl, gr, gl, 'Residual Polarity', rdbu);
plot_row_on_surf(3, 4, Resid_Node_ZCR, sr, sl, gr, gl, 'Residual ZCR', parula(256));
plot_row_on_surf(4, 4, Resid_Node_Poly, sr, sl, gr, gl, 'Residual Poly Degree', cmap_hot);

%% --- HELPER FUNCTIONS ---
function plot_row_on_surf(row_idx, total_rows, plot_vals, sr, sl, gr, gl, tit, cmap)
    cr = zeros(size(gr.cdata)); cl = zeros(size(gl.cdata));
    cr(gr.cdata ~= 0) = plot_vals(gr.cdata(gr.cdata ~= 0));
    cl(gl.cdata ~= 0) = plot_vals(gl.cdata(gl.cdata ~= 0));
    
    clims = [min(plot_vals), max(plot_vals)];
    if diff(clims) == 0, clims = [0 1]; end
    if contains(tit, 'Polarity'), clims = [-1 1]; end 
    
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