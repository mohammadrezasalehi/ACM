% Analyzes how combining sessions improves Test-Retest reliability.
% Compares 1-session, 2-session (Same Day), 2-session (Diff Day), and 4-session data.

clc; clear; close all;

%% 1. Configuration & Subject Mapping
project_root = 'F:\PhD Code\My_PhD_Project';

% Paths
test_dir_node = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features', 'Node_Level');
test_dir_net  = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features', 'EdgeInNetwork_Level');

retest_dir_node = fullfile(project_root, 'results', 'Modulated_FC', 'retest', 'Aggregated_Features', 'Node_Level');
retest_dir_net  = fullfile(project_root, 'results', 'Modulated_FC', 'retest', 'Aggregated_Features', 'EdgeInNetwork_Level');

% =========================================================================
% IMPORTANT: Enter the row indices of the 8 retest subjects as they appear 
% in the original 100-subject array (AllSubjs_...). 
% For example, if they were subjects #5, #12, #45, etc.
% =========================================================================
shared_subj_indices = [8, 14, 19, 26, 31, 55, 57, 59]; % <--- UPDATE THESE NUMBERS!
n_shared = length(shared_subj_indices);

if n_shared ~= 8
    warning('You mentioned 8 subjects, but provided a different number of indices.');
end

%% 2. Load Data
fprintf('Loading Test (100 subjs) and Retest (8 subjs) Data...\n');

% Load Test Data (Extract only the 8 shared subjects)
T_Node = load(fullfile(test_dir_node, 'AllSubjs_Node_Features.mat'));
T_Net  = load(fullfile(test_dir_net, 'AllSubjs_Net_Features.mat'));

% Load Retest Data (All 8 subjects)
R_Node = load(fullfile(retest_dir_node, 'RetestSubjs_Node_Features.mat'));
R_Net  = load(fullfile(retest_dir_net, 'RetestSubjs_Net_Features.mat'));

%% 3. Reliability Calculation Function
% Computes the 4 dose-response conditions for any given feature set
calc_dose = @(T_cell, R_cell) compute_dose_response(T_cell, R_cell, shared_subj_indices, n_shared);

fprintf('Calculating Dose-Response Reliabilities...\n');

% --- Node Level ---
dose_node_var   = calc_dose(T_Node.feat_node_var,   R_Node.feat_node_var);
dose_node_slope = calc_dose(T_Node.feat_node_slope, R_Node.feat_node_slope);
dose_node_zcr   = calc_dose(T_Node.feat_node_zcr,   R_Node.feat_node_zcr);
dose_node_poly  = calc_dose(T_Node.feat_node_poly,  R_Node.feat_node_poly);

% --- Network Level ---
dose_net_var    = calc_dose(T_Net.feat_net_var,     R_Net.feat_net_var);
dose_net_slope  = calc_dose(T_Net.feat_net_slope,   R_Net.feat_net_slope);
dose_net_zcr    = calc_dose(T_Net.feat_net_zcr,     R_Net.feat_net_zcr);
dose_net_poly   = calc_dose(T_Net.feat_net_poly,    R_Net.feat_net_poly);

%% 4. Visualization
bar_labels = {'1 Session\\newline(1200 TP)', '2 Sessions\\newline(Same Day)',...
    '2 Sessions\\newline(Diff Day)', '4 Sessions\\newline(All Data)'};
colors = [0.2 0.6 0.8; 0.9 0.5 0.1; 0.4 0.7 0.4; 0.8 0.2 0.3]; % Distinct colors for bars

% --- Figure 1: Node Level ---
figure('Name', 'Node-Level Dose Response', 'Position', [50, 50, 1200, 800], 'Color', 'w');
sgtitle('Impact of Data Quantity on Node-Level Reliability (Test vs Retest)', 'FontSize', 16, 'FontWeight', 'bold');
plot_dose_bar(1, dose_node_var, 'Variance (Standard Deviation)', bar_labels, colors);
plot_dose_bar(2, dose_node_slope, 'Slope (Absolute Mean)', bar_labels, colors);
plot_dose_bar(3, dose_node_zcr, 'MCR/ZCR (Mean)', bar_labels, colors);
plot_dose_bar(4, dose_node_poly, 'Poly Degree (Non-Linear %)', bar_labels, colors);

% --- Figure 2: Network Level ---
figure('Name', 'Network-Level Dose Response', 'Position', [100, 100, 1200, 800], 'Color', 'w');
sgtitle('Impact of Data Quantity on Edge-in-Network Reliability (Test vs Retest)', 'FontSize', 16, 'FontWeight', 'bold');
plot_dose_bar(1, dose_net_var, 'Variance (Mean Sensitivity)', bar_labels, colors);
plot_dose_bar(2, dose_net_slope, 'Slope (Absolute Mean Sensitivity)', bar_labels, colors);
plot_dose_bar(3, dose_net_zcr, 'MCR/ZCR (Mean Sensitivity)', bar_labels, colors);
plot_dose_bar(4, dose_net_poly, 'Poly Degree (Non-Linear %)', bar_labels, colors);


%% --- HELPER FUNCTIONS ---

function res = compute_dose_response(T_cell, R_cell, idx_map, n_subs)
    % Outputs a 1x4 vector: [Single_Mean, SameDay_Mean, DiffDay_Mean, All_4]
    
    sub_res = zeros(n_subs, 4);
    
    for i = 1:n_subs
        t_idx = idx_map(i);
        r_idx = i; % Retest cell only has 8 rows
        
        % Flatten all 8 sessions (4 Test, 4 Retest)
        T1 = T_cell{t_idx, 1}(:); T2 = T_cell{t_idx, 2}(:);
        T3 = T_cell{t_idx, 3}(:); T4 = T_cell{t_idx, 4}(:);
        
        R1 = R_cell{r_idx, 1}(:); R2 = R_cell{r_idx, 2}(:);
        R3 = R_cell{r_idx, 3}(:); R4 = R_cell{r_idx, 4}(:);
        
        % 1. Single Session (Mean of all 28 pairwise correlations among 8 sessions)
        all_8 = [T1, T2, T3, T4, R1, R2, R3, R4];
        C1 = corr(all_8, 'Rows', 'complete');
        sub_res(i, 1) = mean(C1(tril(true(size(C1)), -1))); % Lower triangle mean
        
        % 2. Same Day (Mean of LR/RL combinations within days)
        T_Day1 = (T1 + T2) / 2; T_Day2 = (T3 + T4) / 2;
        R_Day1 = (R1 + R2) / 2; R_Day2 = (R3 + R4) / 2;
        C2 = corr([T_Day1, T_Day2, R_Day1, R_Day2], 'Rows', 'complete');
        sub_res(i, 2) = mean(C2(tril(true(size(C2)), -1)));
        
        % 3. Diff Day (Mean of LR/LR and RL/RL across days)
        T_LR = (T1 + T3) / 2; T_RL = (T2 + T4) / 2;
        R_LR = (R1 + R3) / 2; R_RL = (R2 + R4) / 2;
        C3 = corr([T_LR, T_RL, R_LR, R_RL], 'Rows', 'complete');
        sub_res(i, 3) = mean(C3(tril(true(size(C3)), -1)));
        
        % 4. All 4 Sessions (Correlation between Mean of Test and Mean of Retest)
        T_All = mean([T1, T2, T3, T4], 2);
        R_All = mean([R1, R2, R3, R4], 2);
        sub_res(i, 4) = corr(T_All, R_All, 'Rows', 'complete');
    end
    
    % Average across the 8 subjects
    res = mean(sub_res, 1, 'omitnan');
end

function plot_dose_bar(idx, vals, tit, labels, colors)
    subplot(2, 2, idx);
    
    % Plot bars individually to assign specific colors
    hold on;
    for b = 1:4
        bar(b, vals(b), 'FaceColor', colors(b, :), 'EdgeColor', 'k', 'LineWidth', 1);
    end
    
    % Formatting
    ylim([0, max(1, max(vals)*1.2)]); % Ensure y-axis goes up to 1 (or slightly higher if needed)
    ylabel('Mean Reliability (Pearson r)', 'FontSize', 10, 'FontWeight', 'bold');
    title(tit, 'FontSize', 12, 'FontWeight', 'bold');
    
    % X-Axis settings
    xticks(1:4);
    % Use sprintf to allow multiline x-labels (\n)
    xticklabels(cellfun(@sprintf, labels, 'UniformOutput', false));
    set(gca, 'TickDir', 'out', 'FontSize', 10);
    grid on;
    
    % Add values on top of bars
    for b = 1:4
        text(b, vals(b) + 0.05, sprintf('%.2f', vals(b)), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
    hold off;
end