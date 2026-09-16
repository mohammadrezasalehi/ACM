%% Suggested Name: plot_within_node_amplitude_stability.m
% Evaluates the Inter-Subject stability of amplitude bounds FOR EACH NODE.
% Plots the distribution of the Interquartile Range (IQR) across subjects,
% demonstrating that identical bins capture identical absolute amplitudes
% for a specific node, regardless of the subject.

clc; clear; close all;

%% 1. Configuration & Data Load
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST1_LR';
clean_pipe = 'No_clean\';

raw_data_dir = fullfile(project_root, 'data', 'inter_parcellated', 'Schaefer200_Kong17', session);

K = 20; 

files = dir(fullfile(raw_data_dir, '*.mat'));
n_subs = min(100, length(files)); 

%% 2. Process Data (Calculate Mean Amplitude per Bin)
fprintf('Calculating Mean Amplitude per Bin...\n');

% Storage: [Subjects x Nodes x K_Bins]
bin_mean_amplitudes = zeros(n_subs, 200, K);

for i = 1:n_subs
    tmp = load(fullfile(files(i).folder, files(i).name));
    ts = tmp.data_on_atlas'; 
    [ntime, nnodes] = size(ts);
    
    for m_node = 1:nnodes
        mod_signal = zscore(ts(:, m_node)); 
        [sorted_signal, ~] = sort(mod_signal);
        bins_indices = round(linspace(1, ntime + 1, K + 1));
        
        for b = 1:K
            idx_start = bins_indices(b);
            idx_end   = bins_indices(b+1) - 1;
            bin_mean_amplitudes(i, m_node, b) = mean(sorted_signal(idx_start:idx_end));
        end
    end
end

%% 3. Calculate Within-Node Inter-Subject Variability (The "Fair" Metric)
fprintf('Calculating Within-Node Inter-Subject Consistency...\n');

% For each Node (200) and each Bin (20), calculate the Interquartile Range (IQR) 
% of the mean amplitudes across the 20 subjects.
% A low IQR means all subjects have almost the exact same amplitude for that bin in that node.

% Storage: [Nodes x K_Bins]
node_iqr_amplitudes = zeros(200, K);

for m_node = 1:nnodes
    for b = 1:K
        subject_vals_for_this_node_bin = bin_mean_amplitudes(:, m_node, b);
        % Calculate IQR (75th percentile - 25th percentile)
        node_iqr_amplitudes(m_node, b) = iqr(subject_vals_for_this_node_bin);
    end
end

%% 4. Visualization
fig = figure('Color', 'w', 'Position', [150, 150, 1000, 600], 'Name', 'Within-Node Amplitude Consistency');

% Draw Boxplots showing the distribution of the 200 IQRs for each Bin
boxplot(node_iqr_amplitudes, 'Colors', 'k', 'Symbol', '.', 'Widths', 0.5);
hold on;

% Add the median line
medians = median(node_iqr_amplitudes, 1);
plot(1:K, medians, 'r-', 'LineWidth', 2, 'DisplayName', 'Median IQR Across Nodes');

% Add a horizontal line to show a "tightness threshold" (e.g., 0.2 Z-score difference)
yline(0.2, 'b--', 'LineWidth', 1.5, 'Label', '0.2 Z-Score Dispersion (Highly Consistent)', ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');

% Aesthetics
grid on; box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1.2);
xlabel('Quantile Bin Index (1 = Lowest, 20 = Highest)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Inter-Subject Dispersion (IQR of Z-Scores)', 'FontSize', 14, 'FontWeight', 'bold');
title('Within-Node Consistency of Amplitude Bins Across Subjects', 'FontSize', 16, 'FontWeight', 'bold');

% Add explanatory text box
annotation('textbox', [0.15 0.75 0.35 0.1], 'String', ...
    {'Low Interquartile Range (IQR) indicates that for', ...
    'any specific node, the amplitude boundaries of a bin', ...
    'are nearly identical across different individuals.'}, ...
    'EdgeColor', 'k', 'FontSize', 11, 'BackgroundColor', 'w');

hold off;