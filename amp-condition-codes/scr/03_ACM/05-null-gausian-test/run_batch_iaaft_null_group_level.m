% Generates 1 IAAFT surrogate per subject for all 100 subjects.
% Aggregates data to calculate Group-Level correlation with Real Data.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR'; 
clean_pipe = 'No_clean\';

raw_data_dir = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', session);
real_group_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Group_Level', session);

% Load Atlas Info for Network Level mapping
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17; n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
net_idx_map = zeros(n_rois, 1);
for i=1:n_rois, net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i})); end

%% 2. Process All Subjects (Generate IAAFT & Extract Features)
files = dir(fullfile(raw_data_dir, '*.mat'));
n_subs = length(files);

% Pre-allocate accumulators for Grand Mean (Surrogate)
Sum_Surr_Net_Var  = zeros(n_nets, n_nets);
Sum_Surr_Net_Slop = zeros(n_nets, n_nets);
Sum_Surr_Net_ZCR  = zeros(n_nets, n_nets);
Sum_Surr_Net_Poly = zeros(n_nets, n_nets);

fprintf('Generating 1 IAAFT Surrogate for %d subjects...\n', n_subs);
h = waitbar(0, 'Processing IAAFT Surrogates...');

for i = 1:n_subs
    % Load Raw Data
    tmp = load(fullfile(files(i).folder, files(i).name));
    ts_real = tmp.data_on_atlas'; 
    [T, N] = size(ts_real);
    
    % Generate 1 Phase-Randomized Surrogate
    F_real = fft(ts_real);
    Amp = abs(F_real);
    Phase = angle(F_real);
    
    rand_phase = 2 * pi * rand(floor(T/2), 1);
    if mod(T,2) == 0
        full_rand_phase = [0; rand_phase; -flipud(rand_phase(1:end-1))];
    else
        full_rand_phase = [0; rand_phase; -flipud(rand_phase)];
    end
    Phase_surr = Phase + repmat(full_rand_phase, 1, N);
    
    F_surr = Amp .* exp(1i * Phase_surr);
    ts_surr = real(ifft(F_surr));
    
    % Run Dynamic Pipeline
    S_out = calc_dynamic_modulation_metrics(ts_surr, 'n_levels', 20);
    
    % Aggregate to Network Level & Add to Sum
    Sum_Surr_Net_Var  = Sum_Surr_Net_Var  + map_edges_to_network(mean(S_out.variance, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    Sum_Surr_Net_Slop = Sum_Surr_Net_Slop + map_edges_to_network(mean(abs(S_out.slope), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    Sum_Surr_Net_ZCR  = Sum_Surr_Net_ZCR  + map_edges_to_network(mean(S_out.zcr, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    
    % Poly
    is_nonlin_S = (S_out.best_poly_degree == 2 | S_out.best_poly_degree == 3);
    valid_mods_S = sum(S_out.best_poly_degree ~= 255, 2);
    p_prev_edge_S = (sum(is_nonlin_S, 2) ./ max(valid_mods_S, 1)) * 100;
    Sum_Surr_Net_Poly = Sum_Surr_Net_Poly + map_edges_to_network(p_prev_edge_S, u_map, v_map, net_idx_map, n_nets);
    
    waitbar(i/n_subs, h);
end
close(h);

% Calculate Grand Mean for Surrogate
G_Surr_Net_Var  = Sum_Surr_Net_Var / n_subs;
G_Surr_Net_Slop = Sum_Surr_Net_Slop / n_subs;
G_Surr_Net_ZCR  = Sum_Surr_Net_ZCR / n_subs;
G_Surr_Net_Poly = Sum_Surr_Net_Poly / n_subs;

%% 3. Load Real Group Data & Calculate Correlations
fprintf('Calculating Group-Level Correlations (Real vs IAAFT)...\n');

% Note: Adjust these filenames if they differ in your Group_Level folder
R_Var  = load(fullfile(real_group_dir, 'Group_Variance_Consensus.mat'));  r_g_var  = R_Var.Consensus_Percent(:);
R_Slop = load(fullfile(real_group_dir, 'Group_Slope_Tmap.mat'));          r_g_slop = R_Slop.T_map(:);
R_ZCR  = load(fullfile(real_group_dir, 'Group_ZCR_Median.mat'));          r_g_zcr  = R_ZCR.Median(:);
R_Poly = load(fullfile(real_group_dir, 'Group_Poly_Prevalence.mat'));     r_g_poly = R_Poly.Prevalence_Percent(:);

% Map Real Group to Network Level (17x17)
R_Net_Var  = map_edges_to_network(mean(reshape(r_g_var, [], n_rois), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
R_Net_Slop = map_edges_to_network(mean(abs(reshape(r_g_slop, [], n_rois)), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
R_Net_ZCR  = map_edges_to_network(mean(reshape(r_g_zcr, [], n_rois), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
% Assume prevalence is already % across subjects, just map it
R_Net_Poly = map_edges_to_network(mean(reshape(r_g_poly, [], n_rois), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);

% Calculate Final Pearson Correlations (Output for your Bar Chart - Row 4, Col 3)
fprintf('\n--- IAAFT GROUP LEVEL CORRELATIONS ---\n');
fprintf('Variance: r = %.4f\n', corr(R_Net_Var(:), G_Surr_Net_Var(:), 'Rows', 'complete'));
fprintf('Slope:    r = %.4f\n', corr(R_Net_Slop(:), G_Surr_Net_Slop(:), 'Rows', 'complete'));
fprintf('ZCR/MCR:  r = %.4f\n', corr(R_Net_ZCR(:), G_Surr_Net_ZCR(:), 'Rows', 'complete'));
fprintf('Poly Deg: r = %.4f\n', corr(R_Net_Poly(:), G_Surr_Net_Poly(:), 'Rows', 'complete'));

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