%% Suggested Name: run_subject_level_residualization.m
% 1. Runs Exact Triplet-Level Spatial Regression (OLS) for each subject.
% 2. Calculates R2 for FC-Only and Full (FC+AR1) models.
% 3. Saves the Pure Dynamic Residuals (19900 x 200) for each subject separately.

clc; clear; close all;

%% 1. Configuration & Paths
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

% Input Directories
dyn_dir  = fullfile(project_root, 'results', 'Modulated_FC', session, clean_pipe);
pred_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Control_Predictors', session, clean_pipe);

% Output Directory (Subject-wise Residuals)
resid_out_dir = fullfile(project_root, 'results', 'Modulated_FC_Residuals', session, clean_pipe);
if ~exist(resid_out_dir, 'dir'), mkdir(resid_out_dir); end

%% 2. Load Predictors (Only once, they contain all subjects)
fprintf('Loading Predictor Matrices (This may take a minute)...\n');
pred_file = fullfile(pred_dir, 'AllSubjs_Static_Predictors_Triplet.mat');
P = load(pred_file);

dyn_files = dir(fullfile(dyn_dir, '*_dyn_mod.mat'));
n_subs = length(dyn_files);

if length(P.FC_XY) ~= n_subs
    error('Mismatch between Predictors and Dynamic subjects!');
end

[n_edges, n_rois] = size(P.FC_XY{1});
n_elements = n_edges * n_rois;

metrics = {'variance', 'slope', 'zcr', 'best_poly_degree'};

% Storage for R^2 values
R2_Summary = struct();
for m = 1:length(metrics)
    R2_Summary.(metrics{m}).R2_FC   = NaN(n_subs, 1);
    R2_Summary.(metrics{m}).R2_Full = NaN(n_subs, 1);
end

%% 3. Subject-wise Processing Loop
fprintf('Starting Triplet-Level Regression for %d subjects...\n', n_subs);

for i = 1:n_subs
    fprintf('Processing Subject %d/%d: %s\n', i, n_subs, dyn_files(i).name);
    
    % 1. Load Dynamic Data
    D = load(fullfile(dyn_dir, dyn_files(i).name));
    
    % 2. Extract Predictors for this subject
    FC_XY = double(P.FC_XY{i}(:));
    FC_XZ = double(P.FC_XZ{i}(:));
    FC_YZ = double(P.FC_YZ{i}(:));
    
    AR_X  = double(P.AR_X{i}(:));
    AR_Y  = double(P.AR_Y{i}(:));
    AR_Z  = double(P.AR_Z{i}(:));
    
    % 3. Standardize Predictors (Crucial for OLS)
    X_FC_norm   = zscore_safe([FC_XY, FC_XZ, FC_YZ]);
    X_Full_norm = zscore_safe([FC_XY, FC_XZ, FC_YZ, AR_X, AR_Y, AR_Z]);
    
    % Add Intercept
    X1 = [ones(n_elements, 1), X_FC_norm];
    X2 = [ones(n_elements, 1), X_Full_norm];
    
    % 4. Pre-allocate Subject Residual Struct
    SubjResid = struct();
    
    % 5. Run Regression for each Metric
    for m = 1:length(metrics)
        met = metrics{m};
        Y = double(D.(met)(:));
        
        % Mask out NaNs and Invalid Poly (255)
        valid = isfinite(Y) & all(isfinite(X1), 2) & all(isfinite(X2), 2);
        if strcmp(met, 'best_poly_degree')
            valid = valid & (Y ~= 255);
        end
        
        X1_v = X1(valid, :);
        X2_v = X2(valid, :);
        Y_v  = Y(valid);
        
        % Calculate Total Sum of Squares (SST)
        SST = sum((Y_v - mean(Y_v)).^2);
        
        % --- Model 1: FC Only ---
        beta_FC = (X1_v' * X1_v) \ (X1_v' * Y_v);
        resid_FC = Y_v - (X1_v * beta_FC);
        SSE_FC = sum(resid_FC.^2);
        
        % --- Model 2: Full (FC + AR) ---
        beta_Full = (X2_v' * X2_v) \ (X2_v' * Y_v);
        resid_Full = Y_v - (X2_v * beta_Full);
        SSE_Full = sum(resid_Full.^2);
        
        % Calculate R^2
        if SST > 0
            R2_Summary.(met).R2_FC(i)   = 1 - (SSE_FC / SST);
            R2_Summary.(met).R2_Full(i) = 1 - (SSE_Full / SST);
        end
        
        % Store FULL Residuals back into 19900x200 matrix
        tmp_Full = NaN(n_elements, 1, 'single');
        tmp_Full(valid) = single(resid_Full);
        
        SubjResid.(met) = reshape(tmp_Full, n_edges, n_rois);
    end
    
    % 6. Save Residuals for this subject
    [~, name_no_ext, ~] = fileparts(dyn_files(i).name);
    save_path = fullfile(resid_out_dir, [name_no_ext '_residualized.mat']);
    save(save_path, '-struct', 'SubjResid');
end

% Save R^2 Summary
save(fullfile(project_root, 'results', 'Modulated_FC_Residuals', 'R2_Summary_Statistics.mat'), '-struct', 'R2_Summary');
fprintf('\nAll Subjects Processed and Saved!\n');

%% --- HELPER FUNCTION ---
function Xz = zscore_safe(X)
    Xz = zeros(size(X));
    for j = 1:size(X,2)
        x = X(:,j);
        mu = mean(x, 'omitnan');
        sd = std(x, 0, 'omitnan');
        if sd == 0 || ~isfinite(sd)
            Xz(:,j) = 0;
        else
            Xz(:,j) = (x - mu) ./ sd;
        end
    end
end