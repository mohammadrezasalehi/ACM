%% run_icc_comprehensive.m
% 1. Calculates ICC for Raw Edges (Memory-Optimized).
% 2. Loads previously computed Node & Network data.
% 3. Plots a comprehensive 3-bar comparison for the manuscript.

clc; clear; close all;

%% 1. Configuration & Paths
project_root = 'F:\PhD Code\My_PhD_Project';
sessions = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
clean_pipe = 'No_clean\';
base_res_dir = fullfile(project_root, 'results', 'Modulated_FC');

% Paths for aggregated data
node_dir = fullfile(base_res_dir, 'Aggregated_Features', 'Node_Level');
net_dir  = fullfile(base_res_dir, 'Aggregated_Features', 'EdgeInNetwork_Level');

% Get Subject List
files = dir(fullfile(base_res_dir, sessions{1}, clean_pipe, '*_dyn_mod.mat'));
subj_names = {files.name};
n_subs = length(subj_names);
n_sess = 4;
n_raw_feats = 19900 * 200; % 3,980,000

% Output Matrix: [4 Metrics x 3 Levels (Raw/Node/Net) x 2 (Mean/CI)]
ICC_Results = zeros(4, 3, 2); 

% Helper function for Mean and CI
calc_mean_ci = @(icc_vec) [mean(icc_vec, 'omitnan'), 1.96 * std(icc_vec, 'omitnan') / sqrt(sum(~isnan(icc_vec)))];

%% 2. Process RAW EDGE LEVEL (Memory-Optimized Loop)
fprintf('--- PART 1: RAW EDGE LEVEL (Memory Optimized) ---\n');
metrics = {'variance', 'slope', 'zcr', 'best_poly_degree'};

for m = 1:4
    metric_name = metrics{m};
    fprintf('Loading and calculating ICC for RAW %s...\n', upper(metric_name));
    
    % Pre-allocate single precision array (~6 GB RAM footprint)
    Raw_Data = zeros(n_raw_feats, n_subs, n_sess, 'single');
    
    for i = 1:n_subs
        for s = 1:n_sess
            file_path = fullfile(base_res_dir, sessions{s}, clean_pipe, subj_names{i});
            tmp = load(file_path, metric_name);
            val = tmp.(metric_name)(:);
            
            % Handle polynomial outliers
            if m == 4
                val(val == 255) = NaN;
            end
            
            Raw_Data(:, i, s) = val;
        end
    end
    
    % Compute Vectorized ICC
    fprintf('  Computing ICC matrix...\n');
    raw_icc_vec = compute_icc_21_vectorized(Raw_Data);
    
    % Store Mean and CI
    ICC_Results(m, 1, :) = calc_mean_ci(raw_icc_vec);
    
    % VERY IMPORTANT: Clear memory before loading the next metric
    clear Raw_Data raw_icc_vec tmp val;
    fprintf('  Done. Memory Cleared.\n');
end

%% 3. Process NODE and NETWORK LEVELS
fprintf('\n--- PART 2: NODE & NETWORK LEVELS ---\n');
NodeData = load(fullfile(node_dir, 'AllSubjs_Node_Features.mat'));
NetData  = load(fullfile(net_dir, 'AllSubjs_Net_Features.mat'));

extract_3d = @(cell_mat, n_feats) reshape(cell2mat(cell_mat'), [n_feats, n_sess, n_subs]);
format_3d = @(cell_mat, n_feats) permute(extract_3d(cell_mat, n_feats), [1, 3, 2]);

% Node
Data_Node_Var = format_3d(NodeData.feat_node_var, 200);
Data_Node_Slp = format_3d(NodeData.feat_node_slope, 200);
Data_Node_ZCR = format_3d(NodeData.feat_node_zcr, 200);
Data_Node_Pol = format_3d(NodeData.feat_node_poly, 200); Data_Node_Pol(Data_Node_Pol==255)=NaN;

% Network
Data_Net_Var = format_3d(NetData.feat_net_var, 289);
Data_Net_Slp = format_3d(NetData.feat_net_slope, 289);
Data_Net_ZCR = format_3d(NetData.feat_net_zcr, 289);
Data_Net_Pol = format_3d(NetData.feat_net_poly, 289); Data_Net_Pol(Data_Net_Pol==255)=NaN;

fprintf('Computing ICCs...\n');
% Save to output matrix (Node = Col 2, Net = Col 3)
ICC_Results(1, 2, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Node_Var));
ICC_Results(2, 2, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Node_Slp));
ICC_Results(3, 2, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Node_ZCR));
ICC_Results(4, 2, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Node_Pol));

ICC_Results(1, 3, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Net_Var));
ICC_Results(2, 3, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Net_Slp));
ICC_Results(3, 3, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Net_ZCR));
ICC_Results(4, 3, :) = calc_mean_ci(compute_icc_21_vectorized(Data_Net_Pol));

%% 4. Visualization: The Final 3-Bar Plot
fprintf('\n--- PART 3: RENDERING PUBLICATION FIGURE ---\n');

titles = {'Variance (Modulatory Sensitivity)', 'Slope (Directional Control)', ...
          'MCR/ZCR (Volatility)', 'Polynomial Degree (Non-Linear)'};
x_labels = {'Raw Edge', 'Node Level', 'Edge-in-Network'};

% Colors matching the paper's aesthetic (Blue, Green, Orange)
colors = [0.1 0.4 0.7;   % Raw Edge (Blue)
          0.2 0.6 0.5;   % Node Level (Green)
          0.8 0.4 0.1];  % Network (Orange)

figure('Name', 'Comprehensive ICC Reliability', 'Position', [50, 50, 1400, 800], 'Color', 'w');
sgtitle('Dose-Response Effect of Spatial Aggregation on Reliability (ICC 2,1)', 'FontSize', 16, 'FontWeight', 'bold');

for m = 1:4
    subplot(2, 2, m);
    hold on; grid on;
    
    means = [ICC_Results(m, 1, 1), ICC_Results(m, 2, 1), ICC_Results(m, 3, 1)];
    cis   = [ICC_Results(m, 1, 2), ICC_Results(m, 2, 2), ICC_Results(m, 3, 2)];
    
    % Draw Bars
    b = bar(1:3, means, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.2, 'BarWidth', 0.6);
    b.CData = colors;
    
    % Draw Error Bars
    errorbar(1:3, means, cis, 'k.', 'LineWidth', 1.5, 'CapSize', 8);
    
    % Text values on top
    for j = 1:3
        text(j, means(j) + cis(j) + 0.03, sprintf('%.2f', means(j)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
    
    % Aesthetics
    title(titles{m}, 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Mean ICC (2,1)', 'FontSize', 12, 'FontWeight', 'bold');
    xticks(1:3);
    xticklabels(x_labels);
    
    % Dynamic Y-limit to keep it clean
    ylim([0, max(0.9, max(means + cis) * 1.3)]);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
end

fprintf('All done! The figure is ready for export.\n');


%% ================= HELPER FUNCTION ================= %%
function icc = compute_icc_21_vectorized(data)
    % COMPUTE_ICC_21_VECTORIZED Computes ICC(2,1) extremely fast.
    
    [v, n, k] = size(data);
    
    mean_tot = mean(mean(data, 3, 'omitnan'), 2, 'omitnan'); 
    mean_sub = mean(data, 3, 'omitnan');                     
    mean_ses = mean(data, 2, 'omitnan');                     
    
    SST  = sum(sum((data - mean_tot).^2, 3, 'omitnan'), 2, 'omitnan'); 
    BSS  = k * sum((mean_sub - mean_tot).^2, 2, 'omitnan');
    WSS  = SST - BSS;
    BSSj = n * sum((mean_ses - mean_tot).^2, 3, 'omitnan');
    ESS  = WSS - BSSj;
    
    BMS = BSS / (n - 1);
    JMS = BSSj / (k - 1);
    EMS = ESS / ((n - 1) * (k - 1));
    
    % ICC(2,1) Formula: Absolute Agreement
    icc = (BMS - EMS) ./ (BMS + (k - 1) .* EMS + (k / n) .* (JMS - EMS));
    icc(isnan(icc)) = 0;
end