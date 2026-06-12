% Calculates Subject-Level correlations (Raw, Node, Network) between 
% Real Data and Random Binning Null Data.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_RL'; 

real_dir = fullfile(project_root, 'results', 'Modulated_FC', session, 'No_clean\');
% UPDATE THIS PATH to where your Random Null individual subjects are saved!
rand_dir = fullfile(project_root, 'results', 'Modulated_FC_RandomNull', session, 'No_clean\');

atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17; n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
net_idx_map = zeros(n_rois, 1);
for i=1:n_rois, net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i})); end

%% 2. Loop Through Subjects
files = dir(fullfile(real_dir, '*_dyn_mod.mat'));
n_subs = length(files);

% Accumulators for Mean Correlation [Subjects x 4 Metrics]
c_raw_var = zeros(n_subs,1); c_raw_slop = zeros(n_subs,1); c_raw_zcr = zeros(n_subs,1); c_raw_poly = zeros(n_subs,1);
c_nod_var = zeros(n_subs,1); c_nod_slop = zeros(n_subs,1); c_nod_zcr = zeros(n_subs,1); c_nod_poly = zeros(n_subs,1);
c_net_var = zeros(n_subs,1); c_net_slop = zeros(n_subs,1); c_net_zcr = zeros(n_subs,1); c_net_poly = zeros(n_subs,1);

fprintf('Calculating Subject-Level Correlations for %d subjects...\n', n_subs);
h = waitbar(0, 'Comparing Real vs Random Null...');

for i = 1:n_subs
    % Load Real & Random Null
    subj_name = files(i).name;
    R = load(fullfile(real_dir, subj_name));
    % Adjust random null filename if necessary
    rand_name = strrep(subj_name, '_dyn_mod.mat', '_random_null.mat');
    N = load(fullfile(rand_dir, rand_name)); 
    
    % --- 1. Raw Level Correlation ---
    idx_v = ~isnan(R.variance(:)) & ~isnan(N.variance(:));
    c_raw_var(i) = corr(R.variance(idx_v), N.variance(idx_v));
    
    idx_t = ~isnan(R.slope(:)) & ~isnan(N.slope(:));
    c_raw_slop(i) = corr(R.slope(idx_t), N.slope(idx_t));
    
    idx_z = ~isnan(R.zcr(:)) & ~isnan(N.zcr(:));
    c_raw_zcr(i) = corr(R.zcr(idx_z), N.zcr(idx_z));
    
    idx_p = (R.best_poly_degree(:) ~= 255) & (N.best_poly_degree(:) ~= 255);
    c_raw_poly(i) = corr(double(R.best_poly_degree(idx_p)), double(N.best_poly_degree(idx_p)));
    
    % --- 2. Node Level Aggregation & Correlation ---
    R_v_nod = std(R.variance, 0, 1, 'omitnan')'; N_v_nod = std(N.variance, 0, 1, 'omitnan')';
    c_nod_var(i) = corr(R_v_nod, N_v_nod, 'Rows', 'complete');
    
    R_t_nod = mean(abs(R.slope), 1, 'omitnan')'; N_t_nod = mean(abs(N.slope), 1, 'omitnan')';
    c_nod_slop(i) = corr(R_t_nod, N_t_nod, 'Rows', 'complete');
    
    R_z_nod = mean(R.zcr, 1, 'omitnan')'; N_z_nod = mean(N.zcr, 1, 'omitnan')';
    c_nod_zcr(i) = corr(R_z_nod, N_z_nod, 'Rows', 'complete');
    
    R_p_nod = (sum(R.best_poly_degree==2 | R.best_poly_degree==3, 1) ./ sum(R.best_poly_degree~=255, 1))';
    N_p_nod = (sum(N.best_poly_degree==2 | N.best_poly_degree==3, 1) ./ sum(N.best_poly_degree~=255, 1))';
    c_nod_poly(i) = corr(R_p_nod, N_p_nod, 'Rows', 'complete');
    
    % --- 3. Network Level Aggregation & Correlation ---
    R_v_net = map_edges_to_network(mean(R.variance, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    N_v_net = map_edges_to_network(mean(N.variance, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    c_net_var(i) = corr(R_v_net(:), N_v_net(:), 'Rows', 'complete');
    
    R_t_net = map_edges_to_network(mean(abs(R.slope), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    N_t_net = map_edges_to_network(mean(abs(N.slope), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    c_net_slop(i) = corr(R_t_net(:), N_t_net(:), 'Rows', 'complete');
    
    R_z_net = map_edges_to_network(mean(R.zcr, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    N_z_net = map_edges_to_network(mean(N.zcr, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    c_net_zcr(i) = corr(R_z_net(:), N_z_net(:), 'Rows', 'complete');
    
    R_p_edge = (sum(R.best_poly_degree==2 | R.best_poly_degree==3, 2) ./ sum(R.best_poly_degree~=255, 2));
    N_p_edge = (sum(N.best_poly_degree==2 | N.best_poly_degree==3, 2) ./ sum(N.best_poly_degree~=255, 2));
    R_p_net = map_edges_to_network(R_p_edge, u_map, v_map, net_idx_map, n_nets);
    N_p_net = map_edges_to_network(N_p_edge, u_map, v_map, net_idx_map, n_nets);
    c_net_poly(i) = corr(R_p_net(:), N_p_net(:), 'Rows', 'complete');
    
    waitbar(i/n_subs, h);
end
close(h);

%% 3. Print Results
fprintf('\n--- RANDOM BINNING SUBJECT-LEVEL CORRELATIONS (Mean over %d Subjs) ---\n', n_subs);
fprintf('RAW LEVEL:     Var=%.4f, Slop=%.4f, ZCR=%.4f, Poly=%.4f\n', mean(c_raw_var), mean(c_raw_slop), mean(c_raw_zcr), mean(c_raw_poly));
fprintf('NODE LEVEL:    Var=%.4f, Slop=%.4f, ZCR=%.4f, Poly=%.4f\n', mean(c_nod_var), mean(c_nod_slop), mean(c_nod_zcr), mean(c_nod_poly));
fprintf('NETWORK LEVEL: Var=%.4f, Slop=%.4f, ZCR=%.4f, Poly=%.4f\n', mean(c_net_var), mean(c_net_slop), mean(c_net_zcr), mean(c_net_poly));

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