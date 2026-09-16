%% run_icc_retest_dose_response.m
% Analyzes how combining sessions improves Test-Retest reliability (ICC 2,1).
% Evaluates: 1-session, 2-session (Same Day), 2-session (Diff Day), and 4-session data.

clc; clear; close all;

%% 1. Configuration & Paths
project_root = 'F:\PhD Code\My_PhD_Project';

test_dir_node = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features', 'Node_Level');
test_dir_net  = fullfile(project_root, 'results', 'Modulated_FC', 'Aggregated_Features', 'EdgeInNetwork_Level');

retest_dir_node = fullfile(project_root, 'results', 'Modulated_FC', 'retest', 'Aggregated_Features', 'Node_Level');
retest_dir_net  = fullfile(project_root, 'results', 'Modulated_FC', 'retest', 'Aggregated_Features', 'EdgeInNetwork_Level');

% Indices of the 8 retest subjects in the 100-subject array
shared_subj_indices = [8, 14, 19, 26, 31, 55, 57, 59]; 
n_shared = length(shared_subj_indices);

%% 2. Load Data
fprintf('Loading Test and Retest Data...\n');
T_Node = load(fullfile(test_dir_node, 'AllSubjs_Node_Features.mat'));
T_Net  = load(fullfile(test_dir_net, 'AllSubjs_Net_Features.mat'));

R_Node = load(fullfile(retest_dir_node, 'RetestSubjs_Node_Features.mat'));
R_Net  = load(fullfile(retest_dir_net, 'RetestSubjs_Net_Features.mat'));

%% 3. Reliability Calculation
fprintf('Calculating Dose-Response ICC(2,1)...\n');

% Calculate for all features
[dose_node_var, ci_node_var]     = compute_dose_icc(T_Node.feat_node_var, R_Node.feat_node_var, shared_subj_indices);
[dose_node_slope, ci_node_slope] = compute_dose_icc(T_Node.feat_node_slope, R_Node.feat_node_slope, shared_subj_indices);
[dose_node_zcr, ci_node_zcr]     = compute_dose_icc(T_Node.feat_node_zcr, R_Node.feat_node_zcr, shared_subj_indices);
[dose_node_poly, ci_node_poly]   = compute_dose_icc(T_Node.feat_node_poly, R_Node.feat_node_poly, shared_subj_indices);

[dose_net_var, ci_net_var]       = compute_dose_icc(T_Net.feat_net_var, R_Net.feat_net_var, shared_subj_indices);
[dose_net_slope, ci_net_slope]   = compute_dose_icc(T_Net.feat_net_slope, R_Net.feat_net_slope, shared_subj_indices);
[dose_net_zcr, ci_net_zcr]       = compute_dose_icc(T_Net.feat_net_zcr, R_Net.feat_net_zcr, shared_subj_indices);
[dose_net_poly, ci_net_poly]     = compute_dose_icc(T_Net.feat_net_poly, R_Net.feat_net_poly, shared_subj_indices);

%% 4. Visualization
bar_labels = {'1 Session\\newline(1200 TP)', '2 Sessions\\newline(Same Day)',...
    '2 Sessions\\newline(Diff Day)', '4 Sessions\\newline(All Data)'};
colors = [0.2 0.6 0.8; 0.9 0.5 0.1; 0.4 0.7 0.4; 0.8 0.2 0.3]; 

% --- Figure 1: Node Level ---
figure('Name', 'Node-Level ICC Dose Response', 'Position', [50, 50, 1200, 800], 'Color', 'w');
sgtitle('Impact of Temporal Integration on Node-Level Reliability (ICC 2,1)', 'FontSize', 16, 'FontWeight', 'bold');
plot_dose_bar(1, dose_node_var, ci_node_var, 'Variance', bar_labels, colors);
ylabel('Absolute Agreement (ICC)', 'FontSize', 11, 'FontWeight', 'bold');
plot_dose_bar(2, dose_node_slope, ci_node_slope, 'Slope', bar_labels, colors);
plot_dose_bar(3, dose_node_zcr, ci_node_zcr, 'ZCR', bar_labels, colors);
plot_dose_bar(4, dose_node_poly, ci_node_poly, 'Poly Degree', bar_labels, colors);

% --- Figure 2: Network Level ---
figure('Name', 'Network-Level ICC Dose Response', 'Position', [100, 100, 1200, 800], 'Color', 'w');
sgtitle('Impact of Temporal Integration on Edge-in-Network Reliability (ICC 2,1)', 'FontSize', 16, 'FontWeight', 'bold');
plot_dose_bar(1, dose_net_var, ci_net_var, 'Variance', bar_labels, colors);
ylabel('Absolute Agreement (ICC)', 'FontSize', 11, 'FontWeight', 'bold');
plot_dose_bar(2, dose_net_slope, ci_net_slope, 'Slope', bar_labels, colors);
plot_dose_bar(3, dose_net_zcr, ci_net_zcr, 'ZCR', bar_labels, colors);
plot_dose_bar(4, dose_net_poly, ci_net_poly, 'Poly Degree', bar_labels, colors);

fprintf('Done!\n');


%% --- HELPER FUNCTIONS ---

function [mean_iccs, ci_iccs] = compute_dose_icc(T_cell, R_cell, idx_map)
    % Outputs a 1x4 vector for ICCs and a 1x4 vector for half-width CIs
    % Scenarios: [Single_Mean, SameDay_Mean, DiffDay_Mean, All_4]
    
    n_subs = length(idx_map);
    n_feats = numel(T_cell{1,1});
    
    % We will flatten the data into massive 2D matrices [ (Subs*Feats) x Sessions ]
    % This treats the entire brain topography as a single multivariate "Fingerprint".
    
    % Pre-allocate full matrices
    T1 = zeros(n_feats * n_subs, 1); T2 = zeros(n_feats * n_subs, 1);
    T3 = zeros(n_feats * n_subs, 1); T4 = zeros(n_feats * n_subs, 1);
    R1 = zeros(n_feats * n_subs, 1); R2 = zeros(n_feats * n_subs, 1);
    R3 = zeros(n_feats * n_subs, 1); R4 = zeros(n_feats * n_subs, 1);
    
    idx = 1;
    for i = 1:n_subs
        t_idx = idx_map(i);
        r_idx = i;
        
        extract = @(cell, r, c) reshape(double(cell{r,c}), [], 1);
        
        t1 = extract(T_cell, t_idx, 1); t2 = extract(T_cell, t_idx, 2);
        t3 = extract(T_cell, t_idx, 3); t4 = extract(T_cell, t_idx, 4);
        
        r1 = extract(R_cell, r_idx, 1); r2 = extract(R_cell, r_idx, 2);
        r3 = extract(R_cell, r_idx, 3); r4 = extract(R_cell, r_idx, 4);
        
        % Replace 255 with NaN for polynomial features
        vars = {t1, t2, t3, t4, r1, r2, r3, r4};
        for v=1:8, vars{v}(vars{v}==255) = NaN; end
        [t1, t2, t3, t4, r1, r2, r3, r4] = vars{:};
        
        len = n_feats;
        T1(idx : idx+len-1) = t1; T2(idx : idx+len-1) = t2;
        T3(idx : idx+len-1) = t3; T4(idx : idx+len-1) = t4;
        
        R1(idx : idx+len-1) = r1; R2(idx : idx+len-1) = r2;
        R3(idx : idx+len-1) = r3; R4(idx : idx+len-1) = r4;
        
        idx = idx + len;
    end
    
    mean_iccs = zeros(1, 4);
    ci_iccs   = zeros(1, 4); % We store the "half-width" of the CI to use in errorbar()
    
    % Condition 1: 8 Single Sessions
    Cond1 = [T1, T2, T3, T4, R1, R2, R3, R4];
    [mean_iccs(1), ci_l, ci_u] = calc_true_icc31(Cond1);
    ci_iccs(1) = (ci_u - ci_l) / 2;
    
    % Condition 2: Same Day (Avg LR+RL within same day)
    Cond2 = [(T1 + T2)/2, (T3 + T4)/2, (R1 + R2)/2, (R3 + R4)/2];
    [mean_iccs(2), ci_l, ci_u] = calc_true_icc31(Cond2);
    ci_iccs(2) = (ci_u - ci_l) / 2;
    
    % Condition 3: Diff Day (Avg LR1+LR2, RL1+RL2)
    Cond3 = [(T1 + T3)/2, (T2 + T4)/2, (R1 + R3)/2, (R2 + R4)/2];
    [mean_iccs(3), ci_l, ci_u] = calc_true_icc31(Cond3);
    ci_iccs(3) = (ci_u - ci_l) / 2;
    
    % Condition 4: Test vs Retest (Avg all 4 sessions)
    Cond4 = [mean([T1, T2, T3, T4], 2, 'omitnan'), mean([R1, R2, R3, R4], 2, 'omitnan')];
    [mean_iccs(4), ci_l, ci_u] = calc_true_icc31(Cond4);
    ci_iccs(4) = (ci_u - ci_l) / 2;
end

function plot_dose_bar(idx, vals, cis, tit, labels, colors)
    subplot(2, 4, idx);
    hold on; grid on;
    
    for b = 1:4
        bar(b, vals(b), 'FaceColor', colors(b, :), 'EdgeColor', 'k', 'LineWidth', 1.2, 'BarWidth', 0.6);
%         errorbar(b, vals(b), cis(b), 'k', 'LineWidth', 1.5, 'CapSize', 8);
        text(b, vals(b) + cis(b) + 0.05, sprintf('%.2f', vals(b)), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
    
    % Formatting
    ylim([0, max(1, max(vals + cis)*1.2)]);
%     ylabel('Absolute Agreement (ICC)', 'FontSize', 11, 'FontWeight', 'bold');
    title(tit, 'FontSize', 13, 'FontWeight', 'bold');
    xticks(1:4);
    xticklabels(cellfun(@sprintf, labels, 'UniformOutput', false));
    xtickangle(90);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 8);
    hold off;
end

function [icc, ci_lower, ci_upper] = calc_true_icc31(X)
    % Calculates ICC(3,1) (Two-way mixed effects, absolute agreement)
    % Excludes rows with any NaNs to ensure balanced ANOVA design
    valid_rows = ~any(isnan(X), 2);
    X = X(valid_rows, :);
    
    [n, k] = size(X);
    
    mean_target = mean(X, 2);
    mean_sess = mean(X, 1);
    mean_grand = mean(X(:));
    
    df_n = n - 1;
    df_k = k - 1;
    df_e = df_n * df_k;
    
    % Sum of Squares
    SST = sum((X(:) - mean_grand).^2);
    SSR = k * sum((mean_target - mean_grand).^2); 
    SSC = n * sum((mean_sess - mean_grand).^2);   
    SSE = SST - SSR - SSC;                        
    
    % Mean Squares
    MSR = SSR / df_n;
    MSE = SSE / df_e;
    
    % Calculate ICC(3,1)
    icc = (MSR - MSE) / (MSR + (k - 1) * MSE);
    if isnan(icc) || icc < 0, icc = 0; end
    
    % Calculate exact 95% Confidence Interval for ICC(3,1)
    alpha = 0.05;
    F_stat = MSR / MSE;
    if isnan(F_stat) || F_stat <= 0, F_stat = 1; end
    
    F_alpha_lower = finv(1 - alpha/2, df_n, df_e);
    F_alpha_upper = finv(alpha/2, df_n, df_e);
    
    L_bound = (F_stat / F_alpha_lower - 1) / ((F_stat / F_alpha_lower) + k - 1);
    U_bound = (F_stat / F_alpha_upper - 1) / ((F_stat / F_alpha_upper) + k - 1);
    
    ci_lower = max(0, L_bound);
    ci_upper = min(1, U_bound);
end