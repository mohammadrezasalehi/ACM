%% =========================================================================
% SCRIPT: analyze_retest_reliability_bars.m
% =========================================================================
% DESCRIPTION:
%   Reliability analysis for HCP RETEST subjects.
%
%   For each feature type:
%       1) Mean pairwise reliability across ALL 8 sessions
%       2) Same-day reliability
%       3) Cross-day reliability
%       4) Test vs Retest reliability
%
%   Produces:
%       - 4 bar plots for NODE level
%       - 4 bar plots for NETWORK level
%
% AUTHOR: Mohammad Reza Salehi
% =========================================================================

clc;
clear;
close all;

%% 1. CONFIGURATION
% =========================================================================

project_root = 'F:\PhD Code\My_PhD_Project';

base_dir = fullfile(project_root, ...
    'results', ...
    'Modulated_FC', ...
    'Aggregated_Features');

% -------------------------------------------------------------------------
% LOAD FEATURES
% -------------------------------------------------------------------------

fprintf('Loading feature files...\n');

NodeData = load(fullfile(base_dir, ...
    'Node_Level', ...
    'AllSubjs_Node_Features.mat'));

NetData = load(fullfile(base_dir, ...
    'EdgeInNetwork_Level', ...
    'AllSubjs_Net_Features.mat'));

%% 2. DEFINE TEST / RETEST STRUCTURE
% =========================================================================
%
% Assumed order:
%
% TEST:
%   1 = REST1_LR
%   2 = REST1_RL
%   3 = REST2_LR
%   4 = REST2_RL
%
% RETEST:
%   5 = REST1_LR
%   6 = REST1_RL
%   7 = REST2_LR
%   8 = REST2_RL
%
% IMPORTANT:
% Your loaded cell arrays must contain 8 sessions per subject.
%
% =========================================================================

n_subs = size(NodeData.feat_node_var, 1);
n_sess = size(NodeData.feat_node_var, 2);

if n_sess ~= 8
    error('Expected 8 sessions per subject (4 test + 4 retest)');
end

%% 3. CALCULATE RELIABILITY SCORES
% =========================================================================

fprintf('Calculating NODE-level reliability...\n');

node_var   = calc_all_metrics(NodeData.feat_node_var);
node_slope = calc_all_metrics(NodeData.feat_node_slope);
node_zcr   = calc_all_metrics(NodeData.feat_node_zcr);
node_poly  = calc_all_metrics(NodeData.feat_node_poly);

fprintf('Calculating NETWORK-level reliability...\n');

net_var   = calc_all_metrics(NetData.feat_net_var);
net_slope = calc_all_metrics(NetData.feat_net_slope);
net_zcr   = calc_all_metrics(NetData.feat_net_zcr);
net_poly  = calc_all_metrics(NetData.feat_net_poly);

%% 4. VISUALIZATION - NODE LEVEL
% =========================================================================

labels = { ...
    'All Pairwise', ...
    'Same Day Avg', ...
    'Cross Day Avg', ...
    'Test vs Retest'};

figure( ...
    'Color', 'w', ...
    'Position', [100 100 1200 900], ...
    'Name', 'Node-Level Retest Reliability');

sgtitle( ...
    'Node-Level Reliability Comparison', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

plot_bar_subplot(1, node_var, labels, ...
    'Variance (STD)');

plot_bar_subplot(2, node_slope, labels, ...
    'Slope (|Mean|)');

plot_bar_subplot(3, node_zcr, labels, ...
    'ZCR/MCR');

plot_bar_subplot(4, node_poly, labels, ...
    'Non-Linear Prevalence');

%% 5. VISUALIZATION - NETWORK LEVEL
% =========================================================================

figure( ...
    'Color', 'w', ...
    'Position', [150 150 1200 900], ...
    'Name', 'Network-Level Retest Reliability');

sgtitle( ...
    'Edge-in-Network Reliability Comparison', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

plot_bar_subplot(1, net_var, labels, ...
    'Variance Sensitivity');

plot_bar_subplot(2, net_slope, labels, ...
    'Slope Sensitivity');

plot_bar_subplot(3, net_zcr, labels, ...
    'ZCR Sensitivity');

plot_bar_subplot(4, net_poly, labels, ...
    'Non-Linear Sensitivity');

fprintf('\nAnalysis Complete!\n');

%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================

function scores = calc_all_metrics(feat_cell)

    n_subs = size(feat_cell, 1);

    % -------------------------------------------------------------
    % SESSION INDICES
    % -------------------------------------------------------------

    test_idx   = 1:4;
    retest_idx = 5:8;

    % Same-day pairs
    same_pairs = [ ...
        1 2;
        3 4;
        5 6;
        7 8];

    % Cross-day pairs
    cross_pairs = [ ...
        1 3;
        2 4;
        5 7;
        6 8];

    % -------------------------------------------------------------
    % OUTPUT STORAGE
    % -------------------------------------------------------------

    all_pair_vals   = zeros(n_subs, 1);
    same_day_vals   = zeros(n_subs, 1);
    cross_day_vals  = zeros(n_subs, 1);
    test_retest_vals = zeros(n_subs, 1);

    % =============================================================
    % SUBJECT LOOP
    % =============================================================

    for i = 1:n_subs

        % ---------------------------------------------------------
        % Load sessions
        % ---------------------------------------------------------

        S = cell(1, 8);

        for s = 1:8
            S{s} = feat_cell{i, s}(:);
        end

        % =========================================================
        % 1) ALL PAIRWISE RELIABILITY
        % =========================================================

        pair_corrs = [];

        for a = 1:8
            for b = a+1:8

                r = corr( ...
                    S{a}, ...
                    S{b}, ...
                    'Rows', 'complete', ...
                    'Type', 'Spearman');

                pair_corrs(end+1) = r;
            end
        end

        all_pair_vals(i) = mean(pair_corrs, 'omitnan');

        % =========================================================
        % 2) SAME-DAY RELIABILITY
        % =========================================================
        %
        % Average LR/RL within same day
        %
        % Day1_Test   = mean(1,2)
        % Day2_Test   = mean(3,4)
        % Day1_Retest = mean(5,6)
        % Day2_Retest = mean(7,8)
        %
        % Then correlate:
        %   Test vs Retest for same day
        %
        % =========================================================

        D1_test   = (S{1} + S{2}) / 2;
        D2_test   = (S{3} + S{4}) / 2;

        D1_retest = (S{5} + S{6}) / 2;
        D2_retest = (S{7} + S{8}) / 2;

        r1 = corr(D1_test, D1_retest, ...
            'Rows', 'complete', ...
            'Type', 'Spearman');

        r2 = corr(D2_test, D2_retest, ...
            'Rows', 'complete', ...
            'Type', 'Spearman');

        same_day_vals(i) = mean([r1 r2], 'omitnan');

        % =========================================================
        % 3) CROSS-DAY RELIABILITY
        % =========================================================
        %
        % Average across days:
        %
        % LR mean = (REST1_LR + REST2_LR)/2
        % RL mean = (REST1_RL + REST2_RL)/2
        %
        % Then correlate LR vs RL
        %
        % =========================================================

        LR_test   = (S{1} + S{3}) / 2;
        RL_test   = (S{2} + S{4}) / 2;

        LR_retest = (S{5} + S{7}) / 2;
        RL_retest = (S{6} + S{8}) / 2;

        r3 = corr(LR_test, RL_test, ...
            'Rows', 'complete', ...
            'Type', 'Spearman');

        r4 = corr(LR_retest, RL_retest, ...
            'Rows', 'complete', ...
            'Type', 'Spearman');

        cross_day_vals(i) = mean([r3 r4], 'omitnan');

        % =========================================================
        % 4) TEST VS RETEST
        % =========================================================

        TEST_mean = mean(cat(2, ...
            S{1}, S{2}, S{3}, S{4}), 2);

        RETEST_mean = mean(cat(2, ...
            S{5}, S{6}, S{7}, S{8}), 2);

        r5 = corr(TEST_mean, RETEST_mean, ...
            'Rows', 'complete', ...
            'Type', 'Spearman');

        test_retest_vals(i) = r5;

    end

    % =============================================================
    % FINAL OUTPUT
    % =============================================================

    scores = [ ...
        mean(all_pair_vals, 'omitnan'), ...
        mean(same_day_vals, 'omitnan'), ...
        mean(cross_day_vals, 'omitnan'), ...
        mean(test_retest_vals, 'omitnan')];

end

%% =========================================================================

function plot_bar_subplot(idx, vals, labels, tit)

    subplot(2,2,idx);

    bar(vals);

    ylim([0 1]);

    xticks(1:4);
    xticklabels(labels);

    xtickangle(20);

    ylabel('Mean Reliability');

    title(tit, ...
        'FontSize', 13, ...
        'FontWeight', 'bold');

    grid on;

    % -------------------------------------------------------------
    % Write numbers on bars
    % -------------------------------------------------------------

    for k = 1:length(vals)

        text( ...
            k, ...
            vals(k) + 0.02, ...
            sprintf('%.3f', vals(k)), ...
            'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', ...
            'FontSize', 11);

    end

end