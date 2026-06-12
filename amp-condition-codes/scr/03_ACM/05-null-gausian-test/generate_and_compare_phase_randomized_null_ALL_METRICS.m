% Generates Phase-Randomized Surrogate Data for Subject 1. 
% Validates the Null model at Raw, Node, and Network levels
% across 100 independent surrogate runs for ALL 4 METRICS.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR'; 
n_runs = 100; % Number of surrogate iterations
p_thresh = 20; % Using 20 levels (bins)

% Paths
raw_data_dir = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', session);
real_res_dir = fullfile(project_root, 'results', 'Modulated_FC', session, 'No_clean\');

% Load Atlas Info
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

% Aggregate Real Results for all 4 metrics
R_raw = struct('v', R.variance, 't', R.slope, 'z', R.zcr, 'p', double(R.best_poly_degree));

R_node = struct(); R_net = struct();
% Node Level Aggregation
R_node.v = std(R_raw.v, 0, 1, 'omitnan')';
R_node.t = mean(abs(R_raw.t), 1, 'omitnan')';
R_node.z = mean(R_raw.z, 1, 'omitnan')';
is_nonlin = (R_raw.p == 2 | R_raw.p == 3);
valid_edges = sum(R.best_poly_degree ~= 255, 1);
R_node.p = (sum(is_nonlin, 1) ./ max(valid_edges, 1))' * 100;

% Network Level Aggregation
R_net.v = map_edges_to_network(mean(R_raw.v, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
R_net.t = map_edges_to_network(mean(abs(R_raw.t), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
R_net.z = map_edges_to_network(mean(R_raw.z, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
valid_mods = sum(R.best_poly_degree ~= 255, 2);
p_prev_edge = (sum(is_nonlin, 2) ./ max(valid_mods, 1)) * 100;
R_net.p = map_edges_to_network(p_prev_edge, u_map, v_map, net_idx_map, n_nets);

%% 3. Generate 100 Surrogate Runs & Run Modulatory Pipeline
fprintf('Generating %d Phase-Randomized Surrogates and computing ALL dynamics...\n', n_runs);

% Pre-allocate arrays to store correlation between Real and Surrogate (4 metrics x 3 levels)
corr_raw  = zeros(n_runs, 4);
corr_node = zeros(n_runs, 4);
corr_net  = zeros(n_runs, 4);

% Prepare FFT of real data
F_real = fft(ts_real);
Amp = abs(F_real);
Phase = angle(F_real);

% Structs to save Run 1 for scatter plots
plot_R = struct('raw', [], 'node', [], 'net', []);
plot_S = struct('raw', [], 'node', [], 'net', []);

h = waitbar(0, 'Processing Surrogate Runs...');
for run = 1:n_runs
    % 1. Create uniform random phase shift
    rand_phase = 2 * pi * rand(floor(T/2), 1);
    if mod(T,2) == 0
        full_rand_phase = [0; rand_phase; -flipud(rand_phase(1:end-1))];
    else
        full_rand_phase = [0; rand_phase; -flipud(rand_phase)];
    end
    Phase_surr = Phase + repmat(full_rand_phase, 1, N);
    
    % 2. Reconstruct surrogate time series
    F_surr = Amp .* exp(1i * Phase_surr);
    ts_surr = real(ifft(F_surr));
    
    % 3. Run full dynamic pipeline
    S_output = calc_dynamic_modulation_metrics(ts_surr, 'n_levels', p_thresh);
    
    % 4. Aggregate Surrogate Results
    S_raw = struct('v', S_output.variance, 't', S_output.slope, 'z', S_output.zcr, 'p', double(S_output.best_poly_degree));
    
    S_node = struct(); S_net = struct();
    S_node.v = std(S_raw.v, 0, 1, 'omitnan')';
    S_node.t = mean(abs(S_raw.t), 1, 'omitnan')';
    S_node.z = mean(S_raw.z, 1, 'omitnan')';
    is_nonlin_S = (S_raw.p == 2 | S_raw.p == 3);
    valid_edges_S = sum(S_output.best_poly_degree ~= 255, 1);
    S_node.p = (sum(is_nonlin_S, 1) ./ max(valid_edges_S, 1))' * 100;
    
    S_net.v = map_edges_to_network(mean(S_raw.v, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    S_net.t = map_edges_to_network(mean(abs(S_raw.t), 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    S_net.z = map_edges_to_network(mean(S_raw.z, 2, 'omitnan'), u_map, v_map, net_idx_map, n_nets);
    valid_mods_S = sum(S_output.best_poly_degree ~= 255, 2);
    p_prev_edge_S = (sum(is_nonlin_S, 2) ./ max(valid_mods_S, 1)) * 100;
    S_net.p = map_edges_to_network(p_prev_edge_S, u_map, v_map, net_idx_map, n_nets);
    
    % 5. Correlate Real vs Surrogate for all 4 metrics
    metrics = {'v', 't', 'z', 'p'};
    for m = 1:4
        field = metrics{m};
        
        idx_raw = ~isnan(R_raw.(field)(:)) & ~isnan(S_raw.(field)(:)) & (R_raw.p(:) ~= 255);
        corr_raw(run, m) = corr(R_raw.(field)(idx_raw), S_raw.(field)(idx_raw));
        
        idx_node = ~isnan(R_node.(field)) & ~isnan(S_node.(field));
        corr_node(run, m) = corr(R_node.(field)(idx_node), S_node.(field)(idx_node));
        
        idx_net = ~isnan(R_net.(field)(:)) & ~isnan(S_net.(field)(:));
        corr_net(run, m) = corr(R_net.(field)(idx_net), S_net.(field)(idx_net));
        
        if run == 1
            plot_R.raw{m} = R_raw.(field)(idx_raw); plot_S.raw{m} = S_raw.(field)(idx_raw);
            plot_R.node{m} = R_node.(field)(idx_node); plot_S.node{m} = S_node.(field)(idx_node);
            plot_R.net{m} = R_net.(field)(idx_net); plot_S.net{m} = S_net.(field)(idx_net);
        end
    end
    
    waitbar(run/n_runs, h);
end
close(h);

%% 4. Visualization (4 Separate Figures)
metric_names = {'Variance', 'Slope (Absolute)', 'MCR / ZCR', 'Polynomial Degree'};

for m = 1:4
    figure('Color', 'w', 'Position', [m*50, m*50, 1200, 800], 'Name', ['Surrogate Validation: ' metric_names{m}]);
    sgtitle(sprintf('Validation against Phase-Randomized Surrogate\nMetric: %s (100 Runs, Subj 1)', metric_names{m}), 'FontSize', 16, 'FontWeight', 'bold');

    % Row 1: Raw
    plot_surrogate_panel(1, plot_R.raw{m}, plot_S.raw{m}, corr_raw(:, m), 'Raw Edge-Level Dynamics');
    
    % Row 2: Node
    plot_surrogate_panel(2, plot_R.node{m}, plot_S.node{m}, corr_node(:, m), 'Node-Level Modulator Profile');
    
    % Row 3: Network
    plot_surrogate_panel(3, plot_R.net{m}, plot_S.net{m}, corr_net(:, m), 'Edge-in-Network Target Profile');
end

%% --- HELPER FUNCTIONS ---
function plot_surrogate_panel(row_idx, R_vec, S_vec, corr_dist, tit)
    % Panel A: Scatter of Run 1
    subplot(3, 2, (row_idx-1)*2 + 1);
    
    % Z-score standardization strictly for visual comparison
    if std(R_vec) > 0, R_norm = (R_vec - mean(R_vec)) / std(R_vec); else, R_norm = R_vec; end
    if std(S_vec) > 0, S_norm = (S_vec - mean(S_vec)) / std(S_vec); else, S_norm = S_vec; end
    
    % FIX 1: Choose plotting method based on data size (row_idx)
    if row_idx == 1
        % Row 1 is Raw data (~4 million points) -> MUST use binscatter
        binscatter(R_norm, S_norm, [150, 150]);
        colormap(gca, 'jet');
    else
        % Rows 2 and 3 are Node (200) and Network (289) -> Use normal scatter
        scatter(R_norm, S_norm, 30, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerFaceColor', [0 0.4470 0.7410]);
    end
    
    hold on; yline(0, 'k--'); xline(0, 'k--');
    
    % FIX 2: Calculate both Pearson and Spearman correlations
    r_pearson  = corr(R_vec, S_vec, 'Type', 'Pearson');
    r_spearman = corr(R_vec, S_vec, 'Type', 'Spearman');
    
    % Display both correlations elegantly on the plot
    txt_str = sprintf('Pearson r = %.3f\nSpearman {\\rho} = %.3f', r_pearson, r_spearman);
    text(double(0.95 * max(R_norm)), double(0.85 * min(S_norm)), txt_str, ...
        'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    
    xlabel('Observed Dynamics (Z)'); ylabel('Surrogate Dynamics (Z)');
    title([tit ' (Single Run)']);
    set(gca, 'TickDir', 'out', 'box', 'off');
    
    % Panel B: Histogram of 100 Runs (Pearson)
    subplot(3, 2, (row_idx-1)*2 + 2);
    histogram(corr_dist, 15, 'FaceColor', [0.8500 0.3250 0.0980], 'EdgeColor', 'k'); % Changed color to differentiate from scatter
    
    xlabel('Pearson Correlation (Observed vs Surrogate)');
    ylabel('Number of Surrogate Datasets');
    title([tit ' (100 Runs Histogram)']);
    
    % Force x-axis limits to show it's clustered around zero
    x_max = max(abs(corr_dist)) * 1.5;
    if x_max == 0 || isnan(x_max), x_max = 0.1; end
    xlim([-x_max, x_max]);
    xline(0, 'k-', 'LineWidth', 2); % Central zero line
    
    % Add mean of the histogram
    mean_corr = mean(corr_dist, 'omitnan');
    xline(mean_corr, 'r--', 'LineWidth', 2);
    legend('Distribution', 'Zero', sprintf('Mean r = %.4f', mean_corr), 'Location', 'best');
    
    set(gca, 'TickDir', 'out', 'box', 'off');
end

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