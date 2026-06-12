%% Suggested Name: map_significant_behavioral_features.m
% Identifies and plots brain nodes and network edges that are significantly
% correlated with the 4 latent behavioral factors.

clc; clear; close all;

%% 1. Configuration & Metric Selection
project_root = 'F:\PhD Code\My_PhD_Project';
base_result_dir = fullfile(project_root, 'results', 'Modulated_FC');
behav_file = fullfile(project_root, 'data', 'bihavioral', 'scores', 'scores_04.csv'); 

% --- SELECT YOUR TARGET DYNAMIC METRIC HERE ---
% Options: 'Variance', 'Slope', 'ZCR', 'PolyDegree'
target_metric = 'PolyDegree'; 

alpha_level = 0.005; % Significance threshold (p < 0.05)

% Load Atlas & Surface Info
surfinfo_path = fullfile(project_root, 'data', 'surf_info', 'surfinfo.mat'); 
load(surfinfo_path); % sr, sl, gr, gl

atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');
addpath(genpath(fullfile(project_root, 'src', '04_utility')));
[~, roi_nets] = load_atlas_info(atlas_info_path);
unique_nets = unique(roi_nets, 'stable'); 
n_nets = 17; n_rois = 200;

%% 2. Load and Align Behavioral Data
fprintf('Aligning Behavioral Data...\n');
files = dir(fullfile(base_result_dir, 'REST1_LR', 'No_clean\', '*_dyn_mod.mat')); % Used just for subject IDs
subject_names = {files.name};
n_subs = length(subject_names);

T_behav = readtable(behav_file);
factor_names = {'WellBeing', 'Cognition', 'Internalizing', 'ProcSpeed', 'SubstanceUse'};
n_factors = length(factor_names);
Y_all = zeros(n_subs, n_factors);

for i = 1:n_subs
    sid_str = regexp(subject_names{i}, '\d{6}', 'match', 'once'); 
    sid = str2double(sid_str);
    idx = find(T_behav.Subject == sid);
    if ~isempty(idx)
        Y_all(i, 1) = T_behav.factor1(idx);
        Y_all(i, 2) = T_behav.factor2(idx);
        Y_all(i, 3) = T_behav.factor3(idx);
        Y_all(i, 4) = T_behav.factor4(idx);
    else
        Y_all(i, :) = NaN;
    end
end

% Remove NaNs
valid_subs = ~isnan(Y_all(:,1));
Y_all = Y_all(valid_subs, :);
n_subs = size(Y_all, 1);

%% 3. Load & Select Dynamic Features
fprintf('Loading Aggregated Features for Metric: %s...\n', target_metric);
NodeData = load(fullfile(base_result_dir, 'Aggregated_Features', 'Node_Level', 'AllSubjs_Node_Features.mat'));
NetData  = load(fullfile(base_result_dir, 'Aggregated_Features', 'EdgeInNetwork_Level', 'AllSubjs_Net_Features.mat'));

% Pre-allocate feature matrices
X_Node = zeros(n_subs, n_rois);
X_Net  = zeros(n_subs, n_nets * n_nets);

sub_orig_idx = find(valid_subs);

% Select the appropriate cell array based on user choice
switch target_metric
    case 'Variance',   C_Node = NodeData.feat_node_var;   C_Net = NetData.feat_net_var;
    case 'Slope',      C_Node = NodeData.feat_node_slope; C_Net = NetData.feat_net_slope;
    case 'ZCR',        C_Node = NodeData.feat_node_zcr;   C_Net = NetData.feat_net_zcr;
    case 'PolyDegree', C_Node = NodeData.feat_node_poly;  C_Net = NetData.feat_net_poly;
end

% Average over the 4 sessions to get Trait-Level features
for i = 1:n_subs
    real_i = sub_orig_idx(i);
    X_Node(i, :) = mean([C_Node{real_i,1}(:), C_Node{real_i,2}(:), C_Node{real_i,3}(:), C_Node{real_i,4}(:)], 2)';
    X_Net(i, :)  = mean([C_Net{real_i,1}(:), C_Net{real_i,2}(:), C_Net{real_i,3}(:), C_Net{real_i,4}(:)], 2)';
end

%% 4. Calculate Correlations & Plot Node-Level (Surface)
fprintf('Calculating Correlations and Rendering Brain Maps...\n');

% Create Custom Blue-White-Red Colormap (Center = 0 = White)
n_c = 256;
c_bwr = [linspace(0,1,n_c/2)', linspace(0,1,n_c/2)', linspace(1,1,n_c/2)'; 
         linspace(1,1,n_c/2)', linspace(1,0,n_c/2)', linspace(1,0,n_c/2)'];

fig_node = figure('Name', 'Significant Nodes', 'Position', [50, 50, 1200, 900], 'Color', 'w');
sgtitle(sprintf('Significant Nodes Associated with Behavior\nMetric: %s (Unthresholded Color, p<%.2f isolated)', target_metric, alpha_level), 'FontSize', 16, 'FontWeight', 'bold');

for f = 1:n_factors
    y_target = Y_all(:, f);
    
    % Correlate each node with behavior
    [r_vec, p_vec] = corr(X_Node, y_target, 'Type', 'Pearson');
    
    % Mask non-significant nodes (set to 0 so they appear white/grey)
    r_vec_masked = r_vec;
    r_vec_masked(p_vec >= alpha_level) = 0;
    
    % Plot on Brain Surface (1 row per behavior)
    plot_row_on_surf(f, n_factors, r_vec_masked', sr, sl, gr, gl, factor_names{f}, c_bwr);
end

%% 5. Calculate Correlations & Plot Network-Level (Heatmap)
fig_net = figure('Name', 'Significant Networks', 'Position', [100, 100, 1000, 800], 'Color', 'w');
sgtitle(sprintf('Significant Edge-in-Network Associations\nMetric: %s (p < %.2f)', target_metric, alpha_level), 'FontSize', 16, 'FontWeight', 'bold');

for f = 1:n_factors
    y_target = Y_all(:, f);
    
    % Correlate each network edge with behavior
    [r_vec, p_vec] = corr(X_Net, y_target, 'Type', 'Pearson');
    
    % Mask non-significant
    r_vec_masked = r_vec;
    r_vec_masked(p_vec >= alpha_level) = 0;
    
    % Reshape back to 17x17
    r_mat_masked = reshape(r_vec_masked, n_nets, n_nets);
    
    % Plot Heatmap
    subplot(2, 2, f);
    imagesc(r_mat_masked);
    colormap(gca, c_bwr); colorbar;
    
    % Force symmetric color scale around 0
    max_val = max(abs(r_mat_masked(:)));
    if max_val == 0, max_val = 0.1; end % Prevent error if nothing is significant
    caxis([-max_val, max_val]);
    
    % Formatting
    title(factor_names{f}, 'FontSize', 14, 'FontWeight', 'bold');
    xticks(1:n_nets); yticks(1:n_nets);
    xticklabels(unique_nets); yticklabels(unique_nets);
    set(gca, 'TickDir', 'out', 'TickLabelInterpreter', 'none', 'FontSize', 8);
    axis square;
end

fprintf('Done!\n');


%% ================= HELPER FUNCTIONS ================= %%
function plot_row_on_surf(row_idx, total_rows, plot_vals, sr, sl, gr, gl, tit, cmap)
    % Map 200 ROIs to vertices
    cr = zeros(size(gr.cdata)); cl = zeros(size(gl.cdata));
    cr(gr.cdata ~= 0) = plot_vals(gr.cdata(gr.cdata ~= 0));
    cl(gl.cdata ~= 0) = plot_vals(gl.cdata(gl.cdata ~= 0));
    
    % Symmetric color limits around 0
    max_val = max(abs(plot_vals));
    if max_val == 0, max_val = 0.1; end
    c_lims = [-max_val, max_val];
    
    is_gifti = isfield(sr, 'data');
    idx_base = (row_idx - 1) * 4;
    
    % View 1: Left Lateral
    ax1 = subplot(total_rows, 4, idx_base + 1);
    plot_surf(sl, cl, is_gifti, -90, c_lims);
    colormap(ax1, cmap); 
    text(ax1, -150, 0, 0, tit, 'FontSize', 14, 'FontWeight', 'bold', ...
         'Interpreter', 'none', 'HorizontalAlignment', 'right');
    
    % View 2: Left Medial
    ax2 = subplot(total_rows, 4, idx_base + 2);
    plot_surf(sl, cl, is_gifti, 90, c_lims);
    colormap(ax2, cmap); 
    
    % View 3: Right Lateral
    ax3 = subplot(total_rows, 4, idx_base + 3);
    plot_surf(sr, cr, is_gifti, 90, c_lims);
    colormap(ax3, cmap); 
    
    % View 4: Right Medial
    ax4 = subplot(total_rows, 4, idx_base + 4);
    plot_surf(sr, cr, is_gifti, -90, c_lims);
    colormap(ax4, cmap); 
    
    % Colorbar
    c = colorbar(ax4, 'Position', [0.92, ax4.Position(2)+0.02, 0.01, ax4.Position(4)-0.04]);
    c.FontSize = 8;
end

function plot_surf(s, c, is_gifti, az, clims)
    if is_gifti
        th = trisurf(s.data{2}.data+1, s.data{1}.data(:,1), s.data{1}.data(:,2), s.data{1}.data(:,3), c);
    else
        th = trisurf(s.faces, s.vertices(:,1), s.vertices(:,2), s.vertices(:,3), c); 
    end
    set(th, 'edgecolor', 'none'); 
    axis image off; 
    set(gca, 'clim', clims);
    view(az, 0); 
    material dull; 
    camlight headlight; 
    lighting gouraud;
end