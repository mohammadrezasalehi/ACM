% Computes and visualizes Test-Retest Reliability across spatial aggregation levels
% (Raw, Node, Network) for both Subject-Level and Group-Level analyses.

clc; clear; close all;

%% 1. Configuration & Paths
project_root = 'F:\PhD Code\My_PhD_Project';
sessions = {'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'};
clean_pipe = 'No_clean\';

base_res_dir = fullfile(project_root, 'results', 'Modulated_FC');
node_dir = fullfile(base_res_dir, 'Aggregated_Features', 'Node_Level');
net_dir  = fullfile(base_res_dir, 'Aggregated_Features', 'EdgeInNetwork_Level');
group_dir = fullfile(base_res_dir, 'Group_Level');

% Load Atlas Info (for Network Level)
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17; n_rois = 200;
[u_map, v_map] = find(triu(true(n_rois), 1));
net_idx_map = zeros(n_rois, 1);
for i=1:n_rois, net_idx_map(i) = find(strcmp(unique_nets, roi_nets{i})); end

% Subject info
files = dir(fullfile(base_res_dir, sessions{1}, clean_pipe, '*_dyn_mod.mat'));
subj_names = {files.name};
n_subs = length(subj_names);
n_sess = length(sessions);

% Output matrix: [4 Metrics x 3 Levels x 2 Types (Subject/Group)]
% Metrics: 1=Var, 2=Slope, 3=ZCR, 4=Poly
Rel_Results = zeros(4, 3, 2);

%% 2. Calculate Reliability: Subject-Level & Group-Level (RAW LEVEL)
fprintf('Calculating RAW LEVEL Reliability (Heavy Computation)...\n');
h = waitbar(0, 'Raw Level Subject Correlations...');
r_raw_subj = zeros(n_subs, 4); % Store mean pairwise corr for each subject

% Accumulators for Group-Level Raw
Sum_Raw_Var = zeros(19900*200, n_sess);
Sum_Raw_Slp = zeros(19900*200, n_sess);
Sum_Raw_ZCR = zeros(19900*200, n_sess);
Sum_Raw_Pol = zeros(19900*200, n_sess);

for i = 1:n_subs
    tmp_v = zeros(19900*200, n_sess); tmp_t = zeros(19900*200, n_sess);
    tmp_z = zeros(19900*200, n_sess); tmp_p = zeros(19900*200, n_sess);
    
    for s = 1:n_sess
        D = load(fullfile(base_res_dir, sessions{s}, clean_pipe, subj_names{i}));
        tmp_v(:,s) = D.variance(:); tmp_t(:,s) = D.slope(:);
        tmp_z(:,s) = D.zcr(:);      tmp_p(:,s) = double(D.best_poly_degree(:));
        
        % Add to group accumulator
        Sum_Raw_Var(:,s) = Sum_Raw_Var(:,s) + tmp_v(:,s);
        Sum_Raw_Slp(:,s) = Sum_Raw_Slp(:,s) + tmp_t(:,s);
        Sum_Raw_ZCR(:,s) = Sum_Raw_ZCR(:,s) + tmp_z(:,s);
        Sum_Raw_Pol(:,s) = Sum_Raw_Pol(:,s) + tmp_p(:,s);
    end
    
    % Subject-Level Mean Pairwise Correlation
    calc_mean_r = @(mat) mean(nonzeros(tril(corr(mat, 'Rows', 'complete'), -1)));
    
    r_raw_subj(i,1) = calc_mean_r(tmp_v);
    r_raw_subj(i,2) = calc_mean_r(tmp_t);
    r_raw_subj(i,3) = calc_mean_r(tmp_z);
    
    % Handle Poly (Ignore 255)
    tmp_p(tmp_p == 255) = NaN;
    r_raw_subj(i,4) = calc_mean_r(tmp_p);
    
    waitbar(i/n_subs, h);
end
close(h);

% 1A. Subject-Level Raw (Mean across 100 subjects)
Rel_Results(:, 1, 1) = mean(r_raw_subj, 1, 'omitnan');

% 1B. Group-Level Raw
calc_mean_r = @(mat) mean(nonzeros(tril(corr(mat, 'Rows', 'complete'), -1)));
Sum_Raw_Pol(Sum_Raw_Pol/n_subs > 10) = NaN; % Naive handling of 255s in sum
Rel_Results(1, 1, 2) = calc_mean_r(Sum_Raw_Var);
Rel_Results(2, 1, 2) = calc_mean_r(Sum_Raw_Slp);
Rel_Results(3, 1, 2) = calc_mean_r(Sum_Raw_ZCR);
Rel_Results(4, 1, 2) = calc_mean_r(Sum_Raw_Pol);

%% 3. Calculate Reliability: Subject-Level (NODE & NETWORK LEVELS)
fprintf('Calculating NODE & NETWORK LEVEL Subject Reliability...\n');
NodeData = load(fullfile(node_dir, 'AllSubjs_Node_Features.mat'));
NetData  = load(fullfile(net_dir, 'AllSubjs_Net_Features.mat'));

r_nod_subj = zeros(n_subs, 4); r_net_subj = zeros(n_subs, 4);

for i = 1:n_subs
    % --- NODE LEVEL ---
    v_n=zeros(200,4); t_n=zeros(200,4); z_n=zeros(200,4); p_n=zeros(200,4);
    % --- NET LEVEL ---
    v_e=zeros(289,4); t_e=zeros(289,4); z_e=zeros(289,4); p_e=zeros(289,4);
    
    for s = 1:4
        v_n(:,s)=NodeData.feat_node_var{i,s}(:); t_n(:,s)=NodeData.feat_node_slope{i,s}(:);
        z_n(:,s)=NodeData.feat_node_zcr{i,s}(:); p_n(:,s)=NodeData.feat_node_poly{i,s}(:);
        
        v_e(:,s)=NetData.feat_net_var{i,s}(:);   t_e(:,s)=NetData.feat_net_slope{i,s}(:);
        z_e(:,s)=NetData.feat_net_zcr{i,s}(:);   p_e(:,s)=NetData.feat_net_poly{i,s}(:);
    end
    
    r_nod_subj(i,:) = [calc_mean_r(v_n), calc_mean_r(t_n), calc_mean_r(z_n), calc_mean_r(p_n)];
    r_net_subj(i,:) = [calc_mean_r(v_e), calc_mean_r(t_e), calc_mean_r(z_e), calc_mean_r(p_e)];
end

% 2A. Subject-Level Node
Rel_Results(:, 2, 1) = mean(r_nod_subj, 1, 'omitnan');
% 3A. Subject-Level Network
Rel_Results(:, 3, 1) = mean(r_net_subj, 1, 'omitnan');

%% 4. Calculate Reliability: Group-Level (NODE & NETWORK LEVELS)
fprintf('Calculating NODE & NETWORK LEVEL Group Reliability...\n');

v_g_n=zeros(200,4); t_g_n=zeros(200,4); z_g_n=zeros(200,4); p_g_n=zeros(200,4);
v_g_e=zeros(289,4); t_g_e=zeros(289,4); z_g_e=zeros(289,4); p_g_e=zeros(289,4);

for s = 1:4
    % Node Level Group
    matrix_2d = NodeData.feat_node_var{:,s}; v_g_n(:,s) = mean(cat(2, matrix_2d(:)), 2);
    matrix_2d = NodeData.feat_node_slope{:,s}; t_g_n(:,s) = mean(cat(2, matrix_2d(:)), 2);
    matrix_2d = NodeData.feat_node_zcr{:,s}; z_g_n(:,s) = mean(cat(2, matrix_2d(:)), 2);
    matrix_2d = NodeData.feat_node_poly{:,s}; p_g_n(:,s) = mean(cat(2, matrix_2d(:)), 2);
    
    % Net Level Group
    matrix_2d = NetData.feat_net_var{:, s}; v_g_e(:,s) = mean(cat(2, matrix_2d(:)), 2);
    matrix_2d = NetData.feat_net_slope{:, s}; t_g_e(:,s) = mean(cat(2, matrix_2d(:)), 2);
    matrix_2d = NetData.feat_net_zcr{:, s}; z_g_e(:,s) = mean(cat(2, matrix_2d(:)), 2);
    matrix_2d = NetData.feat_net_poly{:, s}; p_g_e(:,s) = mean(cat(2, matrix_2d(:)), 2);
end

% 2B. Group-Level Node
Rel_Results(:, 2, 2) = [calc_mean_r(v_g_n), calc_mean_r(t_g_n), calc_mean_r(z_g_n), calc_mean_r(p_g_n)];
% 3B. Group-Level Network
Rel_Results(:, 3, 2) = [calc_mean_r(v_g_e), calc_mean_r(t_g_e), calc_mean_r(z_g_e), calc_mean_r(p_g_e)];


%% 5. Visualization (4-Panel Grouped Bar Charts)
fprintf('Rendering Figures...\n');

titles = {'Variance (Modulatory Sensitivity)', 'Slope (Directional Control)', ...
          'MCR/ZCR (Volatility)', 'Polynomial Degree (Non-Linearity)'};
x_labels = {sprintf('Raw Edge\n(~4M pts)'), sprintf('Node Level\n(200 pts)'),...
    sprintf('Edge-in-Network\n(289 pts)')};
legend_labels = {'Subject-Level (Mean r)', 'Group-Level (Grand Mean r)'};

% Colors: Subject = Blue, Group = Orange
colors = [0.0000 0.4470 0.7410; 
          0.8500 0.3250 0.0980];

fig = figure('Name', 'Spatial Aggregation Reliability', 'Position', [100, 100, 1200, 800], 'Color', 'w');
sgtitle('Impact of Spatial Aggregation on Test-Retest Reliability (4 Sessions)', 'FontSize', 16, 'FontWeight', 'bold');

for m = 1:4
    subplot(2, 2, m);
    
    % Extract data for this metric: [3 Levels x 2 Types]
    plot_data = squeeze(Rel_Results(m, :, :));
    
    b = bar(plot_data, 'grouped', 'EdgeColor', 'k', 'LineWidth', 1);
    b(1).FaceColor = colors(1, :);
    b(2).FaceColor = colors(2, :);
    
    hold on;
    % Add text values on top of bars
    % Add text values on top of the bars for maximum clarity
    ngroups = 3;
    nbars   = 2;
    groupwidth = min(0.8, nbars/(nbars+1.5));
    for col = 1:nbars

        x = (1:ngroups) - groupwidth/2 + ...
            (2*col-1) * groupwidth/(2*nbars);

        y = plot_data(:,col);

        for j = 1:ngroups

            text(x(j), ...
                 y(j) + 0.03*sign(y(j)+eps), ...
                 sprintf('%.2f', y(j)), ...
                 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom', ...
                 'FontSize',9, ...
                 'FontWeight','bold');

        end
    end
    
%     for col = 1:2
%         x_pos = b(col).XData;
%         y_val = b(col).YData;
%         for j = 1:3
%             text(x_pos(j), y_val(j) + 0.05, sprintf('%.2f', y_val(j)), ...
%                 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
%         end
%     end
    
    % Formatting
    title(titles{m}, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Mean Pearson r', 'FontSize', 12, 'FontWeight', 'bold');
    ylim([0, 1.15]); % Leave space for text
    
    % X-Axis settings (Using sprint to enable multiline labels)
    xticklabels({'','','',''});
    yl = ylim;
    for k = 1:3

        text(k, ...
             yl(1)-0.08*(yl(2)-yl(1)), ...
             x_labels{k}, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','top', ...
             'FontSize',11);

    end
%     
    if m == 1
        hLegend = legend(legend_labels, 'Orientation', 'vertical',...
            'Location', 'eastoutside', ...  % ????? ?? southoutside ?? eastoutside
            'FontSize', 11, 'Box', 'on');
        set(hLegend, 'Position', [0.07, 0.9, 0.05, 0.05]);  % [x, y, width, height]
    end
    grid on; hold off;
end