%% Suggested Name: analyze_multilevel_real_vs_gaussian.m
% Validates Real Dynamics against Gaussian Null across 4 distinct levels:
% 1) Subject-Level Raw (19900 x 200)
% 2) Subject-Level Node Aggregated (200 x 1)
% 3) Subject-Level Network Aggregated (17 x 17)
% 4) Group-Level Aggregated

clc; clear; close all;

%% 1. Configuration & Paths
project_root = 'F:\PhD Code\My_PhD_Project';
sessions = {'REST1_LR'};
clean_pipe = 'No_clean\';

% Base Paths
real_base_dir  = fullfile(project_root, 'results', 'Modulated_FC');
gauss_base_dir = fullfile(project_root, 'results', 'Modulated_FC_Gaussian');

% Get Subject List
files = dir(fullfile(real_base_dir, sessions{1}, clean_pipe, '*_dyn_mod.mat'));
subject_names = {files.name};
n_subs = length(subject_names);
n_sess = length(sessions);

%% 2. Pre-allocate Accumulators for Mean Correlations
% Dimensions: [Subjects x Sessions]
r_raw_var   = zeros(n_subs, n_sess); r_raw_slope = zeros(n_subs, n_sess); 
r_raw_zcr   = zeros(n_subs, n_sess); r_raw_poly  = zeros(n_subs, n_sess);

% For visualization, we save the data of Subject 1, Session 1
sample_R_raw = struct(); sample_G_raw = struct();

%% 3. LEVEL 1: Subject-Level Raw (19900 x 200) Analysis
fprintf('Running Level 1: Subject-Level Raw Data Analysis...\n');
h = waitbar(0, 'Level 1: Calculating Raw Correlations (very heavy)...');

for i = 1:n_subs
    for s = 1:n_sess
        % Load Real
        R = load(fullfile(real_base_dir, sessions{s}, clean_pipe, subject_names{i}));
        % Load Gaussian
        G = load(fullfile(gauss_base_dir, sessions{s}, clean_pipe, subject_names{i}));
        
        % Flatten matrices [19900 x 200] -> [3.98 Million x 1]
        v_r = R.variance(:); v_g = G.variance(:);
        t_r = R.slope(:);    t_g = G.slope(:);
        z_r = R.zcr(:);      z_g = G.zcr(:);
        p_r = R.best_poly_degree(:); p_g = G.best_poly_degree(:);
        
        % Filter NaNs
        idx_v = ~isnan(v_r) & ~isnan(v_g);
        idx_t = ~isnan(t_r) & ~isnan(t_g);
        idx_z = ~isnan(z_r) & ~isnan(z_g);
        
        % Treat 255 as NaN for poly
        idx_p = (p_r ~= 255) & (p_g ~= 255);
        
        % Calculate Correlations
        r_raw_var(i,s)   = corr(v_r(idx_v), v_g(idx_v));
        r_raw_slope(i,s) = corr(t_r(idx_t), t_g(idx_t));
        r_raw_zcr(i,s)   = corr(z_r(idx_z), z_g(idx_z));
        r_raw_poly(i,s)  = corr(double(p_r(idx_p)), double(p_g(idx_p)));
        
        % Save Sample for Plotting
        if i == 1 && s == 1
            sample_R_raw.v = v_r(idx_v); sample_G_raw.v = v_g(idx_v);
            sample_R_raw.t = t_r(idx_t); sample_G_raw.t = t_g(idx_t);
            sample_R_raw.z = z_r(idx_z); sample_G_raw.z = z_g(idx_z);
            sample_R_raw.p = double(p_r(idx_p)); sample_G_raw.p = double(p_g(idx_p));
        end
    end
    waitbar(i/n_subs, h);
end
close(h);

% Plot Level 1
mean_r_raw = [mean(r_raw_var(:)), mean(r_raw_slope(:)), mean(r_raw_zcr(:)), mean(r_raw_poly(:))];
plot_multilevel_scatter(1, sample_R_raw, sample_G_raw, mean_r_raw, 'Level 1: Subject-Level Raw (3.98M points)', true);

%% 4. LEVEL 2 & 3: Subject-Level Node and Network Aggregated
fprintf('Running Level 2 & 3: Subject-Level Aggregated Analysis...\n');
R_Node = load(fullfile(real_base_dir, 'Aggregated_Features', 'Node_Level', 'AllSubjs_Node_Features.mat'));
R_Net  = load(fullfile(real_base_dir, 'Aggregated_Features', 'EdgeInNetwork_Level', 'AllSubjs_Net_Features.mat'));

G_Node = load(fullfile(gauss_base_dir, 'Aggregated_Features', 'Node_Level', 'AllSubjs_Node_Features.mat'));
G_Net  = load(fullfile(gauss_base_dir, 'Aggregated_Features', 'EdgeInNetwork_Level', 'AllSubjs_Net_Features.mat'));

% Pre-allocate
r_node = zeros(n_subs, n_sess, 4); % 4 metrics
r_net  = zeros(n_subs, n_sess, 4);

for i = 1:n_subs
    for s = 1:n_sess
        % Node Level (FIX: Added (:) to flatten into column vectors)
        r_node(i,s,1) = corr(R_Node.feat_node_var{i,s}(:),   G_Node.feat_node_var{i,s}(:));
        r_node(i,s,2) = corr(R_Node.feat_node_slope{i,s}(:), G_Node.feat_node_slope{i,s}(:));
        r_node(i,s,3) = corr(R_Node.feat_node_zcr{i,s}(:),   G_Node.feat_node_zcr{i,s}(:));
        r_node(i,s,4) = corr(R_Node.feat_node_poly{i,s}(:),  G_Node.feat_node_poly{i,s}(:));
        
        % Network Level (flatten matrices)
        r_net(i,s,1) = corr(R_Net.feat_net_var{i,s}(:),   G_Net.feat_net_var{i,s}(:));
        r_net(i,s,2) = corr(R_Net.feat_net_slope{i,s}(:), G_Net.feat_net_slope{i,s}(:));
        r_net(i,s,3) = corr(R_Net.feat_net_zcr{i,s}(:),   G_Net.feat_net_zcr{i,s}(:));
        r_net(i,s,4) = corr(R_Net.feat_net_poly{i,s}(:),  G_Net.feat_net_poly{i,s}(:));
    end
end

% Extract Sample 1 for Plotting (FIX: Added (:) here too)
samp_R_node.v = R_Node.feat_node_var{1,1}(:); samp_G_node.v = G_Node.feat_node_var{1,1}(:);
samp_R_node.t = R_Node.feat_node_slope{1,1}(:); samp_G_node.t = G_Node.feat_node_slope{1,1}(:);
samp_R_node.z = R_Node.feat_node_zcr{1,1}(:); samp_G_node.z = G_Node.feat_node_zcr{1,1}(:);
samp_R_node.p = R_Node.feat_node_poly{1,1}(:); samp_G_node.p = G_Node.feat_node_poly{1,1}(:);

samp_R_net.v = R_Net.feat_net_var{1,1}(:); samp_G_net.v = G_Net.feat_net_var{1,1}(:);
samp_R_net.t = R_Net.feat_net_slope{1,1}(:); samp_G_net.t = G_Net.feat_net_slope{1,1}(:);
samp_R_net.z = R_Net.feat_net_zcr{1,1}(:); samp_G_net.z = G_Net.feat_net_zcr{1,1}(:);
samp_R_net.p = R_Net.feat_net_poly{1,1}(:); samp_G_net.p = G_Net.feat_net_poly{1,1}(:);

% Plot Level 2 & 3
mean_r_node = squeeze(mean(mean(r_node, 1), 2));
plot_multilevel_scatter(2, samp_R_node, samp_G_node, mean_r_node, 'Level 2: Subject-Level Node Aggregated (200 pts)', false);

mean_r_net = squeeze(mean(mean(r_net, 1), 2));
plot_multilevel_scatter(3, samp_R_net, samp_G_net, mean_r_net, 'Level 3: Subject-Level Network Aggregated (289 pts)', false);


%% 5. LEVEL 4: Group-Level (Average across all subjects/sessions)
fprintf('Running Level 4: Group-Level Analysis...\n');
get_grand_mean = @(cell_array) mean(cat(3, cell_array{:}), 3, 'omitnan');

% Flatten 100x4 cell array to 1x400 cell array for easy concatenation
flat_R_var = R_Net.feat_net_var(:); flat_G_var = G_Net.feat_net_var(:);
flat_R_slp = R_Net.feat_net_slope(:); flat_G_slp = G_Net.feat_net_slope(:);
flat_R_zcr = R_Net.feat_net_zcr(:); flat_G_zcr = G_Net.feat_net_zcr(:);
flat_R_ply = R_Net.feat_net_poly(:); flat_G_ply = G_Net.feat_net_poly(:);

% Grand Means
grp_R.v = get_grand_mean(flat_R_var); grp_G.v = get_grand_mean(flat_G_var);
grp_R.t = get_grand_mean(flat_R_slp); grp_G.t = get_grand_mean(flat_G_slp);
grp_R.z = get_grand_mean(flat_R_zcr); grp_G.z = get_grand_mean(flat_G_zcr);
grp_R.p = get_grand_mean(flat_R_ply); grp_G.p = get_grand_mean(flat_G_ply);

% Vectorize
grp_R.v = grp_R.v(:); grp_G.v = grp_G.v(:);
grp_R.t = grp_R.t(:); grp_G.t = grp_G.t(:);
grp_R.z = grp_R.z(:); grp_G.z = grp_G.z(:);
grp_R.p = grp_R.p(:); grp_G.p = grp_G.p(:);

grp_r_vals = [corr(grp_R.v, grp_G.v), corr(grp_R.t, grp_G.t), ...
              corr(grp_R.z, grp_G.z), corr(grp_R.p, grp_G.p)];

plot_multilevel_scatter(4, grp_R, grp_G, grp_r_vals, 'Level 4: Group-Level Aggregated (Grand Mean)', false);


%% --- HELPER FUNCTION FOR PLOTTING ---
function plot_multilevel_scatter(fig_num, R_struct, G_struct, r_vals, main_tit, use_binscatter)
    figure('Name', main_tit, 'Position', [fig_num*50, fig_num*50, 1200, 800], 'Color', 'w');
    sgtitle(main_tit, 'FontSize', 16, 'FontWeight', 'bold');
    
    titles = {'Variance', 'Slope (T-val / Abs Mean)', 'MCR/ZCR', 'Poly Degree / Prevalence'};
    fields = {'v', 't', 'z', 'p'};
    
    for k = 1:4
        subplot(2, 2, k);
        x = G_struct.(fields{k});
        y = R_struct.(fields{k});
        
        if use_binscatter
            binscatter(x, y, [150 150]); colormap('jet');
        else
            scatter(x, y, 40, 'filled', 'MarkerFaceAlpha', 0.6);
        end
        
        hold on;
        % Identity Line
        min_all = min([x; y]); max_all = max([x; y]);
        plot([min_all, max_all], [min_all, max_all], 'r--', 'LineWidth', 2);
        
        % Linear Fit
        p_fit = polyfit(x, y, 1);
        x_fit = linspace(min_all, max_all, 100);
        plot(x_fit, polyval(p_fit, x_fit), 'k-', 'LineWidth', 1.5);
        
        title(sprintf('%s\nMean r = %.2f', titles{k}, r_vals(k)), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Gaussian Null (Static FC)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Real Brain Dynamics', 'FontSize', 10, 'FontWeight', 'bold');
        
        if ~use_binscatter, grid on; end
        set(gca, 'TickDir', 'out');
    end
end