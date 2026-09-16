%% Suggested Name: analyze_uniform_binning_nans.m
% Evaluates the percentage of invalid (empty or under-sampled) bins 
% when using the Uniform Amplitude binning strategy.
% Evaluates across different K (bins) and threshold criteria (min frames).
% Ultra-fast: Only relies on indexing, NO correlation calculation.

clc; clear; close all;

%% 1. Configuration & Data Load
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

raw_data_dir = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', session);

% Parameters to test
k_levels = [5, 10, 20];
min_frames_list = 2:10; % Thresholds: minimum frames required to calculate a valid correlation
n_k = length(k_levels);
n_thresh = length(min_frames_list);

% Get Subject List
files = dir(fullfile(raw_data_dir, '*.mat'));
n_subs = min(100, length(files)); % Testing on 20 subjects is statistically sufficient and very fast

%% 2. Process Data (Just Indexing, No FC calculation)
% Storage: [Subjects x Thresholds x K-levels]
invalid_percentages = zeros(n_subs, n_thresh, n_k);

fprintf('Analyzing Uniform Binning vulnerabilities across %d subjects...\n', n_subs);
h = waitbar(0, 'Scanning amplitude distributions...');

for i = 1:n_subs
    % Load Raw BOLD Time Series
    tmp = load(fullfile(files(i).folder, files(i).name));
    ts = tmp.data_on_atlas'; % [1200 x 200]
    [ntime, nnodes] = size(ts);
    
    for k_idx = 1:n_k
        K = k_levels(k_idx);
        
        % We count the total possible bins across all modulators: (200 nodes * K bins)
        total_possible_bins = nnodes * K;
        
        % We will count how many of these bins are invalid for each threshold
        invalid_counts = zeros(n_thresh, 1);
        
        for m_node = 1:nnodes
            mod_signal = ts(:, m_node);
            
            % Uniform Amplitude Edges
            min_val = min(mod_signal);
            max_val = max(mod_signal);
            edges = linspace(min_val, max_val, K + 1);
            edges(1) = -Inf; edges(end) = Inf;
            
            % Discretize
            bin_assignments = discretize(mod_signal, edges);
            
            % Count how many frames fell into each bin
            frames_per_bin = zeros(K, 1);
            for b = 1:K
                frames_per_bin(b) = sum(bin_assignments == b);
            end
            
            % For each threshold, check how many bins failed the criteria
            for t_idx = 1:n_thresh
                thresh = min_frames_list(t_idx);
                invalid_counts(t_idx) = invalid_counts(t_idx) + sum(frames_per_bin < thresh);
            end
        end
        
        % Convert to percentage
        invalid_percentages(i, :, k_idx) = (invalid_counts / total_possible_bins) * 100;
    end
    waitbar(i/n_subs, h);
end
close(h);

% Calculate Mean across subjects
Mean_Invalid_Perc = squeeze(mean(invalid_percentages, 1)); % [n_thresh x n_k]

%% 3. Visualization
fig = figure('Color', 'w', 'Position', [100, 100, 700, 500], 'Name', 'Uniform Binning Failure Rate');

% Colors for the 3 K levels
colors = [0.8500 0.3250 0.0980; % Orange for K=5
          0.4660 0.6740 0.1880; % Green for K=10
          0.0000 0.4470 0.7410];% Blue for K=20

hold on;
plot_handles = zeros(n_k, 1);

for k_idx = 1:n_k
    plot_handles(k_idx) = plot(min_frames_list, Mean_Invalid_Perc(:, k_idx), ...
        '-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(k_idx,:), ...
        'Color', colors(k_idx,:));
end

% Highlight your specific empirical observation (K=20, Thresh > 5)
target_idx = find(min_frames_list == 6); % "greater than 5" means 6 or more
target_val = Mean_Invalid_Perc(target_idx, 3);
plot(6, target_val, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r'); % Red star

text(6, target_val + 2, sprintf('  %.1f%% Invalid\n(K=20, Min>5)', target_val), ...
    'Color', 'r', 'FontWeight', 'bold', 'FontSize', 10);

% Aesthetics
grid on; box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1.2);
xticks(min_frames_list);
xticklabels(arrayfun(@num2str, min_frames_list, 'UniformOutput', false));

xlabel('Minimum Number of Frames Required per Bin', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Percentage of Invalid Bins (%)', 'FontSize', 14, 'FontWeight', 'bold');
title('Vulnerability of Uniform Amplitude Binning', 'FontSize', 16, 'FontWeight', 'bold');

legend(plot_handles, {'K = 5 Bins', 'K = 10 Bins', 'K = 20 Bins'}, 'Location', 'NorthWest', 'FontSize', 12);

hold off;