% Calculates and plots within-subject test-retest reliability across 4 sessions
% for 4 dynamic modulation metrics (Variance, Slope, ZCR, Poly Degree).

clc; clear; close all;

%% 1. Configuration
base_dir = 'F:\PhD Code\My_PhD_Project\results\Modulated_FC\';
sessions = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
clean_pipe = 'No_clean\';

metrics_names = {'variance', 'slope', 'zcr', 'best_poly_degree'};
titles = {'Variance Reliability', 'Slope Reliability', ...
          'Zero-Crossing Rate (ZCR) Reliability', 'Polynomial Degree Reliability'};

n_sess = length(sessions);
n_metrics = length(metrics_names);

%% 2. Find Common Subjects Across All 4 Sessions
fprintf('Finding common subjects across all sessions...\n');
% Start with subjects in the first session
files = dir(fullfile(base_dir, sessions{1}, clean_pipe, '*_dyn_mod.mat'));
common_subs = {files.name};

% Intersect with the other 3 sessions to ensure we only use complete data
for s = 2:n_sess
    s_files = dir(fullfile(base_dir, sessions{s}, clean_pipe, '*_dyn_mod.mat'));
    s_names = {s_files.name};
    common_subs = intersect(common_subs, s_names);
end

n_subs = length(common_subs);
fprintf('Found %d complete subjects.\n', n_subs);

if n_subs == 0
    error('No common subjects found across all 4 sessions. Check paths!');
end

%% 3. Calculate Subject-Level Reliability
% Pre-allocate a 4D array: [session x session x metric x subject]
rel_matrices = zeros(n_sess, n_sess, n_metrics, n_subs);

h = waitbar(0, 'Calculating Reliability...');
for i = 1:n_subs
    % Load all 4 sessions for subject i
    subj_data = cell(n_sess, 1);
    for s = 1:n_sess
        file_path = fullfile(base_dir, sessions{s}, clean_pipe, common_subs{i});
        subj_data{s} = load(file_path);
    end
    
    % Compute correlation for each metric
    for m = 1:n_metrics
        m_name = metrics_names{m};
        
        % Extract metric from all sessions and vectorize
        % Size of tmp_mat: [Number_of_Edges, Number_of_Sessions]
        num_elements = numel(subj_data{1}.(m_name));
        tmp_mat = zeros(num_elements, n_sess);
        
        for s = 1:n_sess
            val = subj_data{s}.(m_name)(:);
            tmp_mat(:, s) = double(val); % Convert to double (important for uint8)
        end
        
        % Find rows that have NaN or 255 (trivial edge flag) in ANY session
        bad_rows = any(isnan(tmp_mat) | tmp_mat == 255, 2);
        
        % Remove invalid edges to compute clean correlation
        clean_mat = tmp_mat(~bad_rows, :);
        
        % Calculate 4x4 correlation matrix for this subject & metric
        if ~isempty(clean_mat)
            rel_matrices(:, :, m, i) = corr(clean_mat);
        else
            rel_matrices(:, :, m, i) = NaN; % Fallback if all data is invalid
        end
    end
    
    waitbar(i/n_subs, h, sprintf('Processing Subject %d/%d', i, n_subs));
end
close(h);

%% 4. Average Across Subjects
% Mean across the 4th dimension (subjects)
mean_rel = mean(rel_matrices, 4, 'omitnan');

%% 5. Visualization (2x2 Subplot Heatmaps)
figure('Name', 'Dynamic Modulation Reliability', 'Position', [100, 100, 1000, 800], 'Color', 'w');

% Labels for X and Y axes (using TeX formatting to prevent underscore issues)
sess_labels = {'R1\_LR', 'R1\_RL', 'R2\_LR', 'R2\_RL'};

for m = 1:n_metrics
    subplot(2, 2, m);
    
    % Plot Heatmap
    imagesc(mean_rel(:, :, m));
    colormap('jet'); % Using 'jet' to match your requested image style
    colorbar;
    
    % Set axis limits dynamically based on data, max is always 1
    caxis([min(mean_rel(:,:,m), [], 'all') * 0.95, 1]); 
    
    % Formatting Axes
    set(gca, 'XTick', 1:n_sess, 'XTickLabel', sess_labels, ...
             'YTick', 1:n_sess, 'YTickLabel', sess_labels, ...
             'FontSize', 10, 'TickDir', 'out', 'Box', 'off');
    
    title(titles{m}, 'FontSize', 12, 'FontWeight', 'bold');
    
    % Overlay Text Values
    for r = 1:n_sess
        for c = 1:n_sess
            val = mean_rel(r, c, m);
            
            % Smart Text Color: White for dark background (blue/dark red), Black for light (yellow/cyan)
            % This is tailored for the 'jet' colormap
            c_limits = caxis;
            normalized_val = (val - c_limits(1)) / (c_limits(2) - c_limits(1));
            if normalized_val < 0.25 || normalized_val > 0.8
                t_color = 'w';
            else
                t_color = 'k';
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
sgtitle('Mean Subject-Level Reliability of Dynamic Modulation Metrics', 'FontSize', 16, 'FontWeight', 'bold');