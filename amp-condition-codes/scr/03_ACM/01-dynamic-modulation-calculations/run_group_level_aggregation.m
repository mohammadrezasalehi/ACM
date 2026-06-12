%% Suggested Name: run_group_level_aggregation.m
% Aggregates subject-level dynamic modulation metrics into group-level maps.
% Applies advanced statistical methods (T-tests, Non-linear Prevalence, Consensus Maps).

clc; clear; close all;

%% 1. Configuration
% Set paths (adjust session name as needed)
session = 'REST1_RL'; 
clean_pipe = 'No_clean\';
base_dir = 'F:\PhD Code\My_PhD_Project\results\Modulated_FC\';
input_dir = fullfile(base_dir, session, clean_pipe);

% Create Group-Level output directory
group_out_dir = fullfile(base_dir, 'Group_Level', session);
if ~exist(group_out_dir, 'dir')
    mkdir(group_out_dir);
end

% Get Subject List
files = dir(fullfile(input_dir, '*_dyn_mod.mat'));
n_subs = length(files);
fprintf('Found %d subjects in %s. Starting Group Aggregation...\n', n_subs, session);

if n_subs == 0
    error('No subjects found. Check your directories.');
end

%% 2. Pre-allocate 3D Data Matrices [Edges x Nodes x Subjects]
% Assuming Schaefer 200 (19900 edges x 200 nodes)
tmp = load(fullfile(input_dir, files(1).name));
[nedges, nnodes] = size(tmp.variance);

var_3d   = zeros(nedges, nnodes, n_subs, 'single');
slope_3d = zeros(nedges, nnodes, n_subs, 'single');
zcr_3d   = zeros(nedges, nnodes, n_subs, 'single');
poly_3d  = zeros(nedges, nnodes, n_subs, 'uint8');

%% 3. Load Data into RAM
h = waitbar(0, 'Loading Subject Data into Memory...');
for i = 1:n_subs
    data = load(fullfile(input_dir, files(i).name));
    var_3d(:, :, i)   = data.variance;
    slope_3d(:, :, i) = data.slope;
    zcr_3d(:, :, i)   = data.zcr;
    poly_3d(:, :, i)  = data.best_poly_degree;
    
    waitbar(i/n_subs, h, sprintf('Loading Subject %d/%d', i, n_subs));
end
close(h);

%% 4. Group-Level Statistics

% ---------------------------------------------------------
% A. SLOPE: 1-Sample T-Test vs. 0 (Vectorized for extreme speed)
% ---------------------------------------------------------
fprintf('Calculating T-Map for Slope...\n');
mu = mean(slope_3d, 3, 'omitnan');
sigma = std(slope_3d, 0, 3, 'omitnan');
n_valid = sum(~isnan(slope_3d), 3);

% Avoid division by zero for edges with no data
sigma(sigma == 0) = NaN; 
t_val = mu ./ (sigma ./ sqrt(n_valid));

% Calculate two-tailed P-value (Requires Statistics Toolbox)
p_val = 2 * (1 - tcdf(abs(t_val), n_valid - 1));

slope_group = struct('T_map', t_val, 'P_map', p_val, 'Mean', mu);
save(fullfile(group_out_dir, 'Group_Slope_Tmap.mat'), '-struct', 'slope_group');

% ---------------------------------------------------------
% B. POLYNOMIAL DEGREE: Mode and Non-Linearity Prevalence
% ---------------------------------------------------------
fprintf('Calculating Polynomial Maps...\n');
% 1. Prevalence Map: % of subjects where degree is 2 or 3
is_non_linear = (poly_3d == 2) | (poly_3d == 3);
valid_poly_counts = sum(poly_3d ~= 255, 3);
prevalence_map = (sum(is_non_linear, 3) ./ valid_poly_counts) * 100;

% 2. Mode (Most frequent degree)
% Convert to single temporarily to use NaN for mode calculation
poly_single = single(poly_3d);
poly_single(poly_single == 255) = NaN;
mode_poly = mode(poly_single, 3);

poly_group = struct('Prevalence_Percent', prevalence_map, 'Mode', mode_poly);
save(fullfile(group_out_dir, 'Group_Poly_Prevalence.mat'), '-struct', 'poly_group');

% ---------------------------------------------------------
% C. VARIANCE: Median and Consensus of Hubs
% ---------------------------------------------------------
fprintf('Calculating Variance Maps...\n');
% 1. Median
median_var = median(var_3d, 3, 'omitnan');

% 2. Consensus Map (Top 10% Modulators per subject)
consensus_counts = zeros(nedges, nnodes, 'single');
valid_var_counts = sum(~isnan(var_3d), 3);

for i = 1:n_subs
    v = var_3d(:, :, i);
    % Find the 90th percentile threshold for this subject (ignoring NaNs)
    p90 = prctile(v(:), 90); 
    
    % Add 1 to edges that are above the 90th percentile
    is_top10 = (v >= p90);
    consensus_counts = consensus_counts + single(is_top10);
end
% Convert count to probability (percentage of subjects where this edge was a top 10% hub)
consensus_map = (consensus_counts ./ valid_var_counts) * 100;

var_group = struct('Median', median_var, 'Consensus_Percent', consensus_map);
save(fullfile(group_out_dir, 'Group_Variance_Consensus.mat'), '-struct', 'var_group');

% ---------------------------------------------------------
% D. ZCR/MCR: Median
% ---------------------------------------------------------
fprintf('Calculating ZCR/MCR Maps...\n');
median_zcr = median(zcr_3d, 3, 'omitnan');

zcr_group = struct('Median', median_zcr);
save(fullfile(group_out_dir, 'Group_ZCR_Median.mat'), '-struct', 'zcr_group');

fprintf('\nSuccess! All Group-level matrices saved in:\n%s\n', group_out_dir);