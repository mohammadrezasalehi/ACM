% Calculates and plots Group-Level test-retest reliability across 4 sessions.
% Correlates the aggregated group maps (Consensus, T-maps, Medians, Prevalence).

clc; clear; close all;

%% 1. Configuration
base_dir = 'F:\PhD Code\My_PhD_Project\results\Modulated_FC\Group_Level\';
sessions = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
n_sess = length(sessions);

% Define the files and fields we want to extract for the 4 subplots
metrics_info = {
    'Group_Variance_Consensus.mat', 'Consensus_Percent',  'Variance (Consensus Map)';
    'Group_Slope_Tmap.mat',         'T_map',              'Slope (T-Value Map)';
    'Group_ZCR_Median.mat',         'Median',             'ZCR/MCR (Median Map)';
    'Group_Poly_Prevalence.mat',    'Prevalence_Percent', 'Poly Degree (Non-Linear Prevalence)'
};

n_metrics = size(metrics_info, 1);

%% 2. Calculate Group-Level Reliability
% Pre-allocate a 3D array: [session x session x metric]
group_rel_matrices = zeros(n_sess, n_sess, n_metrics);

fprintf('Calculating Group-Level Reliability across %d sessions...\n', n_sess);

for m = 1:n_metrics
    file_name = metrics_info{m, 1};
    field_name = metrics_info{m, 2};
    metric_title = metrics_info{m, 3};
    
    % We don't know the exact size yet, so we load the first session to preallocate
    tmp = load(fullfile(base_dir, sessions{1}, file_name));
    num_elements = numel(tmp.(field_name));
    
    % Matrix to hold flattened data for all 4 sessions [Pixels x 4]
    data_all_sess = zeros(num_elements, n_sess, 'double');
    
    % Load data for all 4 sessions
    for s = 1:n_sess
        file_path = fullfile(base_dir, sessions{s}, file_name);
        data = load(file_path);
        
        % Flatten the 2D matrix (Edges x Nodes) into a 1D column
        val = data.(field_name)(:);
        data_all_sess(:, s) = double(val);
    end
    
    % Find rows that have NaN in ANY session (e.g., trivial edges)
    bad_rows = any(isnan(data_all_sess), 2);
    
    % Remove invalid edges to compute clean correlation
    clean_mat = data_all_sess(~bad_rows, :);
    
    % Calculate 4x4 correlation matrix for this group metric
    group_rel_matrices(:, :, m) = corr(clean_mat);
    
    fprintf('Processed: %s\n', metric_title);
end

%% 3. Visualization (2x2 Subplot Heatmaps)
figure('Name', 'Group-Level Dynamic Modulation Reliability', ...
       'Position', [100, 100, 1000, 800], 'Color', 'w');

% Labels for X and Y axes (using TeX formatting to prevent underscore issues)
sess_labels = {'R1\_LR', 'R1\_RL', 'R2\_LR', 'R2\_RL'};

for m = 1:n_metrics
    subplot(2, 2, m);
    
    % Plot Heatmap
    imagesc(group_rel_matrices(:, :, m));
    colormap('jet'); 
    colorbar;
    
    % Dynamic caxis: Group reliability is usually very high. 
    % We set the min to slightly below the lowest value to show contrast.
    min_val = min(group_rel_matrices(:,:,m), [], 'all');
    caxis([min(0.5, min_val * 0.95), 1]); % Ensures lower bound is at least 0.5 if data is high
    
    % Formatting Axes
    set(gca, 'XTick', 1:n_sess, 'XTickLabel', sess_labels, ...
             'YTick', 1:n_sess, 'YTickLabel', sess_labels, ...
             'FontSize', 10, 'TickDir', 'out', 'Box', 'off');
    
    title(metrics_info{m, 3}, 'FontSize', 12, 'FontWeight', 'bold');
    
    % Overlay Text Values
    for r = 1:n_sess
        for c = 1:n_sess
            val = group_rel_matrices(r, c, m);
            
            % Smart Text Color based on colormap intensity
            c_limits = caxis;
            normalized_val = (val - c_limits(1)) / (c_limits(2) - c_limits(1));
            
            if normalized_val < 0.25 || normalized_val > 0.8
                t_color = 'w'; % White text for dark blue / dark red
            else
                t_color = 'k'; % Black text for cyan / yellow
            end
            
            text(c, r, sprintf('%.3f', val), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', t_color, ...
                'FontWeight', 'bold', ...
                'FontSize', 11);
        end
    end
end

% Adjust overall layout
sgtitle('Group-Level Reliability of Dynamic Modulation Metrics (HCP)', ...
        'FontSize', 16, 'FontWeight', 'bold');