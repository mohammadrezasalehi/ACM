%% Suggested Name: extract_static_control_predictors_v2.m
% Extracts Static FC (Triplets) and Autocorrelation (AR1) for all subjects.
% Note: Regional Variance is excluded because input data is Z-scored (Var=1).

clc; clear; close all;

%% 1. Configuration & Paths
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

raw_data_dir = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', session);
out_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Control_Predictors', session, clean_pipe);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

files = dir(fullfile(raw_data_dir, '*.mat'));
n_subs = length(files);
n_rois = 200;

[u_map, v_map] = find(triu(true(n_rois), 1));
n_edges = length(u_map);

%% 2. Process All Subjects
fprintf('Extracting Predictors (FC Triplets and AR1) for %d subjects...\n', n_subs);

Predictors = struct();
Predictors.FC_XY = cell(n_subs, 1);
Predictors.FC_XZ = cell(n_subs, 1);
Predictors.FC_YZ = cell(n_subs, 1);

Predictors.AR_X  = cell(n_subs, 1);
Predictors.AR_Y  = cell(n_subs, 1);
Predictors.AR_Z  = cell(n_subs, 1);

h = waitbar(0, 'Extracting Static FC and AR1...');

for i = 1:n_subs
    % Load Raw BOLD Time Series (Z-scored)
    tmp = load(fullfile(files(i).folder, files(i).name));
    ts_z = tmp.data_on_atlas'; % [1200 x 200]
    
    % --- 1. Autocorrelation (AR1) ---
    reg_ar1 = zeros(n_rois, 1);
    for node = 1:n_rois
        c = corrcoef(ts_z(1:end-1, node), ts_z(2:end, node));
        reg_ar1(node) = c(1,2);
    end
    
    % --- 2. Static FC ---
    static_fc_mat = corr(ts_z);
    
    % --- 3. Construct Predictor Matrices [19900 x 200] ---
    mat_FC_XY = zeros(n_edges, n_rois, 'single');
    mat_FC_XZ = zeros(n_edges, n_rois, 'single');
    mat_FC_YZ = zeros(n_edges, n_rois, 'single');
    
    mat_AR_X  = zeros(n_edges, n_rois, 'single');
    mat_AR_Y  = zeros(n_edges, n_rois, 'single');
    mat_AR_Z  = zeros(n_edges, n_rois, 'single');
    
    static_vec = static_fc_mat(u_map + (v_map-1)*n_rois);
    base_FC_XY = repmat(static_vec, 1, n_rois);
    
    base_AR_X  = repmat(reg_ar1(u_map), 1, n_rois);
    base_AR_Y  = repmat(reg_ar1(v_map), 1, n_rois);
    
    for mod_node = 1:n_rois
        % FC Triplets
        mat_FC_XY(:, mod_node) = base_FC_XY(:, mod_node);
        mat_FC_XZ(:, mod_node) = static_fc_mat(u_map, mod_node);
        mat_FC_YZ(:, mod_node) = static_fc_mat(v_map, mod_node);
        
        % AR1
        mat_AR_X(:, mod_node)  = base_AR_X(:, mod_node);
        mat_AR_Y(:, mod_node)  = base_AR_Y(:, mod_node);
        mat_AR_Z(:, mod_node)  = reg_ar1(mod_node);
    end
    
    % Store
    Predictors.FC_XY{i} = mat_FC_XY;
    Predictors.FC_XZ{i} = mat_FC_XZ;
    Predictors.FC_YZ{i} = mat_FC_YZ;
    
    Predictors.AR_X{i}  = mat_AR_X;
    Predictors.AR_Y{i}  = mat_AR_Y;
    Predictors.AR_Z{i}  = mat_AR_Z;
    
    waitbar(i/n_subs, h);
end
close(h);

%% 3. Save Predictors
fprintf('\nSaving extracted predictors to disk...\n');
save_path = fullfile(out_dir, 'AllSubjs_Static_Predictors_Triplet.mat');
save(save_path, '-struct', 'Predictors', '-v7.3');
fprintf('Successfully saved predictors to:\n%s\n', save_path);