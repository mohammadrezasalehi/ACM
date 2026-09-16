%% Suggested Name: plot_uniform_bin_occupancy.m
% Visualizes the extreme imbalance in frame occupancy across the 20 bins 
% when using the Uniform Amplitude strategy.
% Demonstrates why Quantile binning (fixed 60 TRs per bin) is superior.

clc; clear; close all;

%% 1. Configuration & Data Load
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

raw_data_dir = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', session);

K = 20; % Number of bins

% Get Subject List
files = dir(fullfile(raw_data_dir, '*.mat'));
n_subs = min(20, length(files)); % 20 subjects are enough for this statistical proof

%% 2. Process Data (Calculate Occupancy per Bin)
fprintf('Calculating bin occupancy for Uniform strategy across %d subjects...\n', n_subs);

% Storage: [Subjects x Nodes x K_Bins] -> Will be flattened for plotting
occupancy_data = zeros(n_subs, 200, K);

h = waitbar(0, 'Scanning amplitude distributions...');
for i = 1:n_subs
    tmp = load(fullfile(files(i).folder, files(i).name));
    ts = tmp.data_on_atlas'; % [1200 x 200]
    [ntime, nnodes] = size(ts);
    
    for m_node = 1:nnodes
        mod_signal = ts(:, m_node);
        
        % Uniform Amplitude Edges
        min_val = min(mod_signal);
        max_val = max(mod_signal);
        edges = linspace(min_val, max_val, K + 1);
        edges(1) = -Inf; edges(end) = Inf;
        
        % Discretize
        bin_assignments = discretize(mod_signal, edges);
        
        % Count frames per bin
        for b = 1:K
            occupancy_data(i, m_node, b) = sum(bin_assignments == b);
        end
    end
    waitbar(i/n_subs, h);
end
close(h);

%% 3. Visualization
% Flatten the data: For each of the 20 bins, we have (20 Subs * 200 Nodes) = 4000 observations
flat_occupancy = reshape(occupancy_data, [], K); % [4000 x 20]

fig = figure('Color', 'w', 'Position', [100, 100, 1000, 600], 'Name', 'Uniform Bin Occupancy');

% Draw Boxplots
boxplot(flat_occupancy, 'Colors', 'k', 'Symbol', 'r.');
hold on;

% Calculate and plot the Quantile (Ideal) baseline for comparison
ideal_frames = 1200 / K; % 60 TRs
yline(ideal_frames, 'b--', 'LineWidth', 2, 'Label', sprintf('Quantile Strategy (%d TRs/Bin)', ideal_frames), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'b');

% Highlight the dangerous zone (less than 6 frames)
fill([0, K+1, K+1, 0], [0, 0, 5, 5], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
yline(5, 'r-', 'LineWidth', 1.5, 'Label', 'Danger Zone (<= 5 TRs)', ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'FontSize', 10, 'Color', 'r');

% Aesthetics
grid on; box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1.2);
xlabel('Bin Index (From Lowest to Highest Amplitude)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Number of Frames (TRs) in Bin', 'FontSize', 14, 'FontWeight', 'bold');
title('The Imbalance of Uniform Amplitude Binning (K=20)', 'FontSize', 16, 'FontWeight', 'bold');

% Force Y-axis to show the full extent
ylim([0, max(flat_occupancy(:)) * 1.1]);

hold off;