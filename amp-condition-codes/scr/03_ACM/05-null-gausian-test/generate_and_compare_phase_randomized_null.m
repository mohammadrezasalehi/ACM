%% Suggested Name: generate_and_compare_phase_randomized_null.m
% Generates Phase-Randomized Surrogate Data (similar to Edge-Centric paper)
% for Subject 1. Validates the Null model at Raw, Node, and Network levels
% across 100 independent surrogate runs.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR'; 
n_runs = 100; % Number of surrogate iterations
p_thresh = 20; % Using 20 levels (bins) for your method

% Paths
raw_data_dir = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', session);
real_res_dir = fullfile(project_root, 'results', 'Modulated_FC', session, 'No_clean\');

% Load Atlas Info for Network Level mapping
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17; n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
net_idx_map = zeros(n_rois, 1);
for i=1:n_rois, net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i})); end

%% 2. Load Real Data (Subject 1)
files = dir(fullfile(raw_data_dir, '*.mat'));
subj_name = files(1).name;
tmp = load(fullfile(raw_data_dir, subj_name));
ts_real = tmp.data_on_atlas'; % [1200 Timepoints x 200 Nodes]
[T, N] = size(ts_real);

% Load Real Dynamic Results for Comparison
R = load(fullfile(real_res_dir, [strrep(subj_name, '.mat', ''), '_dyn_mod.mat']));

% Aggregate Real Results (We focus on Variance as the primary metric for brevity)
R_raw = R.variance;
R_node = std(R_raw, 0, 1, 'omitnan')'; % 200x1
R_net  = map_edges_to_network(mean(R_raw, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);

%% 3. Generate 100 Surrogate Runs & Run Modulatory Pipeline
fprintf('Generating %d Phase-Randomized Surrogates and computing dynamics...\n', n_runs);

% Pre-allocate arrays to store correlation between Real and Surrogate
corr_raw  = zeros(n_runs, 1);
corr_node = zeros(n_runs, 1);
corr_net  = zeros(n_runs, 1);

% Prepare FFT of real data
% According to the paper: "random phase between [0; 2pi] added to each frequency bin uniformly across all brain regions"
F_real = fft(ts_real);
Amp = abs(F_real);
Phase = angle(F_real);

h = waitbar(0, 'Processing Surrogate Runs...');
for run = 1:n_runs
    % 1. Create uniform random phase shift for each frequency bin
    % T/2 phases are generated, then mirrored to ensure real-valued IFFT
    rand_phase = 2 * pi * rand(floor(T/2), 1);
    if mod(T,2) == 0
        full_rand_phase = [0; rand_phase; -flipud(rand_phase(1:end-1))];
    else
        full_rand_phase = [0; rand_phase; -flipud(rand_phase)];
    end
    
    % Apply the SAME random phase shift to ALL brain regions (preserves cross-correlation / Static FC)
    Phase_surr = Phase + repmat(full_rand_phase, 1, N);
    
    % 2. Reconstruct surrogate time series
    F_surr = Amp .* exp(1i * Phase_surr);
    ts_surr = real(ifft(F_surr));
    
    % 3. Run your dynamic pipeline on the surrogate data
    % (Assuming calc_dynamic_modulation_metrics is in your path and optimized)
    S = calc_dynamic_modulation_metrics(ts_surr, 'n_levels', p_thresh, 'req_metrics', {'variance'});
    
    % 4. Aggregate Surrogate Results
    S_raw = S.variance;
    S_node = std(S_raw, 0, 1, 'omitnan')';
    S_net  = map_edges_to_network(mean(S_raw, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    
    % 5. Correlate Real vs Surrogate
    % Flatten matrices and filter NaNs
    idx_raw = ~isnan(R_raw(:)) & ~isnan(S_raw(:));
    corr_raw(run) = corr(R_raw(idx_raw), S_raw(idx_raw));
    
    idx_node = ~isnan(R_node) & ~isnan(S_node);
    corr_node(run) = corr(R_node(idx_node), S_node(idx_node));
    
    idx_net = ~isnan(R_net(:)) & ~isnan(S_net(:));
    corr_net(run) = corr(R_net(idx_net), S_net(idx_net));
    
    % Save data from the FIRST run to plot the 2D Scatter (like panel a in paper)
    if run == 1
        plot_R_raw = R_raw(idx_raw); plot_S_raw = S_raw(idx_raw);
        plot_R_node = R_node(idx_node); plot_S_node = S_node(idx_node);
        plot_R_net = R_net(idx_net); plot_S_net = S_net(idx_net);
    end
    
    waitbar(run/n_runs, h);
end
close(h);

%% 4. Visualization (Mimicking the Paper's Format)
figure('Color', 'w', 'Position', [100, 100, 1200, 800], 'Name', 'Surrogate Validation');
sgtitle('Validation against Phase-Randomized Surrogate (100 Runs, Subject 1)', 'FontSize', 16, 'FontWeight', 'bold');

% Plot Row 1: Raw Edge-Level (19900 x 200 points)
plot_surrogate_panel(1, plot_R_raw, plot_S_raw, corr_raw, 'Raw Edge-Level Dynamics');

% Plot Row 2: Node-Level Aggregated (200 points)
plot_surrogate_panel(2, plot_R_node, plot_S_node, corr_node, 'Node-Level Modulator Profile');

% Plot Row 3: Network-Level Aggregated (289 points)
plot_surrogate_panel(3, plot_R_net, plot_S_net, corr_net, 'Edge-in-Network Target Profile');


%% --- HELPER FUNCTIONS ---
function plot_surrogate_panel(row_idx, R_vec, S_vec, corr_dist, tit)
    % Panel A: 2D Binscatter of Run 1
    subplot(3, 2, (row_idx-1)*2 + 1);
    
    % Standardize for visual comparison (optional, but makes axes comparable)
    R_norm = (R_vec - mean(R_vec)) / std(R_vec);
    S_norm = (S_vec - mean(S_vec)) / std(S_vec);
    
    binscatter(R_norm, S_norm, [100, 100]);
    colormap(gca, 'jet');
    hold on; yline(0, 'k--'); xline(0, 'k--');
    
    r_val = corr(R_vec, S_vec);
    text(double(0.9*max(R_norm)), double(0.9*min(S_norm)), sprintf('r = %.5f', r_val), ...
        'HorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold');
    
    xlabel('Observed Dynamics (Z)');
    ylabel('Surrogate Dynamics (Z)');
    title([tit ' (Single Run)']);
    set(gca, 'TickDir', 'out', 'box', 'off');
    
    % Panel B: Histogram of 100 Runs
    subplot(3, 2, (row_idx-1)*2 + 2);
    histogram(corr_dist, 15, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'k');
    
    xlabel('Correlation (Observed vs Surrogate)');
    ylabel('Number of Surrogate Datasets');
    title([tit ' (100 Runs Histogram)']);
    
    % To mimic the paper, force x-axis limits to show it's clustered around zero
    x_max = max(abs(corr_dist)) * 1.5;
    if x_max == 0, x_max = 0.1; end
    xlim([-x_max, x_max]);
    xline(0, 'r--', 'LineWidth', 1.5);
    set(gca, 'TickDir', 'out', 'box', 'off');
end

function net_mat = map_edges_to_network(edge_vec, u, v, net_idx_map, n_nets)
    sum_scores = zeros(n_nets, n_nets); count_scores = zeros(n_nets, n_nets);
    for e = 1:length(u)
        ni = u(e); nj = v(e);
        net_i = net_idx_map(ni); net_j = net_idx_map(nj);
        val = edge_vec(e);
        sum_scores(net_i, net_j) = sum_scores(net_i, net_j) + val;
        count_scores(net_i, net_j) = count_scores(net_i, net_j) + 1;
        if net_i ~= net_j
            sum_scores(net_j, net_i) = sum_scores(net_j, net_i) + val;
            count_scores(net_j, net_i) = count_scores(net_j, net_i) + 1;
        end
    end
    net_mat = sum_scores ./ count_scores;
end