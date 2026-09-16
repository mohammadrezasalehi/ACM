%% Suggested Name: run_spatial_regression_and_residual_maps_v3.m
% 1. Uses Sufficient Statistics (X'X, X'Y) to fit exact OLS on ~4M points.
% 2. Calculates exact R^2 for FC and Delta R^2 for AR1.
% 3. SAFELY aggregates Residuals (using omitnan) to Node and Network levels.

clc; clear; close all;

%% 1. Configuration & Load Data
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

pred_dir  = fullfile(project_root, 'results', 'Modulated_FC', 'Control_Predictors', session);
dyn_dir   = fullfile(project_root, 'results', 'Modulated_FC', session, clean_pipe);

% Load Atlas Info for Network Level mapping
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17; n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
net_idx_map = zeros(n_rois, 1);
for i=1:n_rois, net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i})); end
n_edges = length(u_map);
n_elements = n_edges * n_rois;

fprintf('Loading Predictor Matrices...\n');
PredData = load(fullfile(pred_dir, 'AllSubjs_Static_Predictors_Triplet.mat'));
n_subs = length(PredData.FC_XY);
dyn_files = dir(fullfile(dyn_dir, '*_dyn_mod.mat'));

%% 2. Process Subjects (Nested Models & Exact OLS)
fprintf('Running Exact Spatial Regression...\n');

R2_FC_Values   = zeros(n_subs, 4);
R2_Full_Values = zeros(n_subs, 4);

% Pre-allocate 3D tensors for SAFE aggregation [Elements x Subjects]
All_Resid_Var  = NaN(n_elements, n_subs, 'single');
All_Resid_Slop = NaN(n_elements, n_subs, 'single');
All_Resid_ZCR  = NaN(n_elements, n_subs, 'single');
All_Resid_Poly = NaN(n_elements, n_subs, 'single');

h = waitbar(0, 'Regressing out Static features...');
for i = 1:n_subs
    % 1. Load Dynamic Data (Y)
    D = load(fullfile(dyn_files(i).folder, dyn_files(i).name));
    
    Y_mats = {D.variance(:), D.slope(:), D.zcr(:), single(D.best_poly_degree(:))};
    valid_mask = ~isnan(Y_mats{1}) & ~isnan(Y_mats{2}) & ~isnan(Y_mats{3}) & (Y_mats{4} ~= 255);
    
    % 2. Build Design Matrices (X)
    X_FC = [ones(sum(valid_mask), 1, 'single'), ...
            zscore(PredData.FC_XY{i}(valid_mask)), ...
            zscore(PredData.FC_XZ{i}(valid_mask)), ...
            zscore(PredData.FC_YZ{i}(valid_mask))]; 
            
    X_Full = [X_FC, ...
              zscore(PredData.AR_X{i}(valid_mask)), ...
              zscore(PredData.AR_Y{i}(valid_mask)), ...
              zscore(PredData.AR_Z{i}(valid_mask))]; 
    
    % 3. Fit Models using Sufficient Statistics
    resid_mats = cell(1, 4);
    
    for m = 1:4
        Y_tgt = Y_mats{m}(valid_mask);
        SST = sum((Y_tgt - mean(Y_tgt)).^2);
        
        % Model 1: FC Only
        Beta_FC = (X_FC' * X_FC) \ (X_FC' * Y_tgt);
        SS_res_FC = sum((Y_tgt - X_FC * Beta_FC).^2);
        R2_FC_Values(i, m) = 1 - (SS_res_FC / SST);
        
        % Model 2: Full
        Beta_Full = (X_Full' * X_Full) \ (X_Full' * Y_tgt);
        Resid_Full = Y_tgt - X_Full * Beta_Full;
        SS_res_Full = sum(Resid_Full.^2);
        R2_Full_Values(i, m) = 1 - (SS_res_Full / SST);
        
        % Store Residuals safely for this subject
        Full_Resid = NaN(n_elements, 1, 'single');
        Full_Resid(valid_mask) = Resid_Full;
        resid_mats{m} = Full_Resid;
    end
    
    All_Resid_Var(:, i)  = resid_mats{1};
    All_Resid_Slop(:, i) = resid_mats{2};
    All_Resid_ZCR(:, i)  = resid_mats{3};
    All_Resid_Poly(:, i) = resid_mats{4};
    
    waitbar(i/n_subs, h);
end
close(h);

%% 3. SAFE Spatial Aggregation (Addressing the NaN Bug)
fprintf('Calculating Grand Mean Residuals and Aggregating...\n');

% Group Level (Mean across 100 subjects, safely omitting NaNs)
Grp_Resid_Var  = reshape(mean(All_Resid_Var, 2, 'omitnan'), n_edges, n_rois);
Grp_Resid_Slop = reshape(mean(All_Resid_Slop, 2, 'omitnan'), n_edges, n_rois);
Grp_Resid_ZCR  = reshape(mean(All_Resid_ZCR, 2, 'omitnan'), n_edges, n_rois);
Grp_Resid_Poly = reshape(mean(All_Resid_Poly, 2, 'omitnan'), n_edges, n_rois);

% Node Level Aggregation
Node_Resid_Var  = std(Grp_Resid_Var, 0, 1, 'omitnan')';
Node_Resid_Slop = mean(abs(Grp_Resid_Slop), 1, 'omitnan')';
Node_Resid_ZCR  = mean(Grp_Resid_ZCR, 1, 'omitnan')';
Node_Resid_Poly = mean(Grp_Resid_Poly, 1, 'omitnan')'; 

% Network Level Aggregation (17x17)
Net_Resid_Var  = map_edges_to_network(mean(Grp_Resid_Var, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
Net_Resid_Slop = map_edges_to_network(mean(abs(Grp_Resid_Slop), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
Net_Resid_ZCR  = map_edges_to_network(mean(Grp_Resid_ZCR, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
Net_Resid_Poly = map_edges_to_network(mean(Grp_Resid_Poly, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);

%% 4. Print R^2 Summary
fprintf('\n--- R^2 SUMMARY (Variance Explained by Static FC & AR1) ---\n');
mean_R2_FC = mean(R2_FC_Values, 1);
mean_R2_AR = mean(R2_Full_Values, 1) - mean_R2_FC; % Delta R^2
metrics = {'Variance', 'Slope', 'MCR/ZCR', 'Poly Degree'};
for m = 1:4
    fprintf('%s:\t StaticFC = %.1f%% \t +AR1 (Delta) = %.1f%% \t Total = %.1f%%\n', ...
        metrics{m}, mean_R2_FC(m)*100, mean_R2_AR(m)*100, (mean_R2_FC(m)+mean_R2_AR(m))*100);
end

%% --- HELPER FUNCTION ---
function net_mat = map_edges_to_network(edge_vec, u, v, net_idx_map, n_nets)
    sum_scores = zeros(n_nets, n_nets); count_scores = zeros(n_nets, n_nets);
    for e = 1:length(u)
        ni = u(e); nj = v(e);
        net_i = net_idx_map(ni); net_j = net_idx_map(nj);
        val = edge_vec(e);
        if ~isnan(val)
            sum_scores(net_i, net_j) = sum_scores(net_i, net_j) + val;
            count_scores(net_i, net_j) = count_scores(net_i, net_j) + 1;
            if net_i ~= net_j
                sum_scores(net_j, net_i) = sum_scores(net_j, net_i) + val;
                count_scores(net_j, net_i) = count_scores(net_j, net_i) + 1;
            end
        end
    end
    net_mat = sum_scores ./ count_scores;
end