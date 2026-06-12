% Calculates node-level modulatory features (Mean, Entropy, Polarity, etc.)
% Maps the 200-dimensional profiles onto the 3D cortical surface.

clc; clear; close all;

%% 1. Configuration
project_root = 'F:\PhD Code\My_PhD_Project';
session = 'REST2_LR'; 
group_dir = fullfile(project_root, 'results', 'Modulated_FC', 'Group_Level', session);

surfinfo_path = fullfile(project_root, 'data', 'surf_info', 'surfinfo.mat'); 
atlas_info_path = fullfile(project_root, 'data', 'atlases', 'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');

% Load Surf Info
fprintf('Loading Cortical Surfaces...\n');
load(surfinfo_path); % sr, sl, gr, gl

%% 2. Load Group-Level Matrices [19900 edges x 200 nodes]
fprintf('Loading Group-Level matrices...\n');
V = load(fullfile(group_dir, 'Group_Variance_Consensus.mat')); var_mat = V.Consensus_Percent;
T = load(fullfile(group_dir, 'Group_Slope_Tmap.mat'));         t_mat = T.T_map;
Z = load(fullfile(group_dir, 'Group_ZCR_Median.mat'));         zcr_mat = Z.Median;
P = load(fullfile(group_dir, 'Group_Poly_Prevalence.mat'));    poly_prev = P.Prevalence_Percent;
poly_mode = P.Mode;

[n_edges, n_rois] = size(var_mat);

%% 3. Calculate Node-Level Features (Collapsing Edges)
fprintf('Calculating Node-Level Features...\n');
NodeFeatures = struct();

% --- A. Variance Features ---
var_mat_clean = var_mat; var_mat_clean(isnan(var_mat_clean)) = 0;
NodeFeatures.Var_Mean = mean(var_mat_clean, 1)'; % Global Modulatory Strength
NodeFeatures.Var_HubDeg = sum(var_mat_clean > 50, 1)'; % Count of edges with >50% consensus

p_var = var_mat_clean ./ sum(var_mat_clean, 1);
p_var(p_var == 0) = 1; % log(1)=0, safe for entropy
NodeFeatures.Var_Entropy = -sum(p_var .* log2(p_var), 1)'; % Modulatory Specificity

% --- B. Slope (T-Map) Features ---
NodeFeatures.Slope_Facil = sum(t_mat > 2.5, 1)'; % Facilitator Index
NodeFeatures.Slope_Inhib = sum(t_mat < -2.5, 1)'; % Inhibitor Index
% Polarity: +1 (Pure Facil), -1 (Pure Inhib), 0 (Balanced)
NodeFeatures.Slope_Polarity = (NodeFeatures.Slope_Facil - NodeFeatures.Slope_Inhib) ./ ...
                              (NodeFeatures.Slope_Facil + NodeFeatures.Slope_Inhib + eps);

% --- C. ZCR/MCR Features ---
NodeFeatures.ZCR_Mean = mean(zcr_mat, 1, 'omitnan')'; % Global Destabilizer
zcr_sorted = sort(zcr_mat, 1, 'descend', 'MissingPlacement', 'last');
NodeFeatures.ZCR_Top5 = mean(zcr_sorted(1:5, :), 1, 'omitnan')'; % Top-5 Edge Impact

% --- D. Polynomial (Non-Linear) Features ---
NodeFeatures.Poly_Footprint = mean(poly_prev, 1, 'omitnan')'; % Mean Non-Linearity
% Entropy of Degrees (Complexity)
p0 = sum(poly_mode == 0, 1); p1 = sum(poly_mode == 1, 1);
p2 = sum(poly_mode == 2, 1); p3 = sum(poly_mode == 3, 1);
P_deg = [p0; p1; p2; p3] ./ sum([p0; p1; p2; p3], 1);
P_deg(P_deg == 0) = 1;
NodeFeatures.Poly_Entropy = -sum(P_deg .* log2(P_deg), 1)'; 

% Save all 200x1 vectors to a mat file
save(fullfile(group_dir, 'Node_Level_Features.mat'), '-struct', 'NodeFeatures');
fprintf('Node features saved to Group_Level directory.\n');

%% 4. Visualization on 3D Cortical Surface
% We will create 1 figure per variable, containing multiple rows (one for each metric)

% --- Figure 1: Variance ---
fig1 = figure('Name', 'Node-Level Variance Maps', 'Position', [50, 50, 1200, 900], 'Color', 'w');
sgtitle('Variance Modulatory Features (Consensus Map)', 'FontSize', 18, 'FontWeight', 'bold');
plot_row_on_surf(1, 3, NodeFeatures.Var_Mean, sr, sl, gr, gl,  {'Global Modulatory', 'Strength (Mean)'}, turbo);
plot_row_on_surf(2, 3, NodeFeatures.Var_Entropy, sr, sl, gr, gl, {'Modulatory Diversity', '(Shannon Entropy)'}, turbo);
plot_row_on_surf(3, 3, NodeFeatures.Var_HubDeg, sr, sl, gr, gl, {'Consensus Hub', 'Degree (>50%)'}, hot);

% --- Figure 2: Slope (Facilitation vs Inhibition) [LOG SCALE VERSION] ---
% Custom Red-Blue colormap for Polarity
rdbu = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];

fig2 = figure('Name', 'Node-Level Slope Maps', 'Position', [100, 100, 1200, 900], 'Color', 'w');
sgtitle('Slope / Directional Modulatory Features (Log10 Scale)', 'FontSize', 18, 'FontWeight', 'bold');

% Apply Log10(x+1) transformation for visualization
log_facil = log10(NodeFeatures.Slope_Facil + 1);
log_inhib = log10(NodeFeatures.Slope_Inhib + 1);

% Find global max of the log-transformed data to set a unified colorbar limit
global_log_max = max([max(log_facil), max(log_inhib)]);

% Temporarily override the plot_row_on_surf function's c_lims behavior
% by passing the data and we will rely on a slightly modified logic inside the function 
% or just let the function use the new log data. The function automatically takes [min max].
% To force them to be identical, we will explicitly set the limits AFTER calling the function.

plot_row_on_surf(1, 3, log_facil, sr, sl, gr, gl, {'Facilitator Index', '(Log10 Count)'}, autumn);
plot_row_on_surf(2, 3, log_inhib, sr, sl, gr, gl, {'Inhibitor Index', '(Log10 Count)'}, winter);
plot_row_on_surf(3, 3, NodeFeatures.Slope_Polarity, sr, sl, gr, gl, {'Modulatory Polarity', '(-1=Inhib, +1=Facil)'}, rdbu);

% FORCE UNIFIED COLOR LIMITS FOR ROW 1 AND ROW 2
axes_objs = findobj(fig2, 'Type', 'axes');
for ax = axes_objs'
    % Check if the title of the axis contains 'Log10 Count'
    if contains(ax.Title.String, 'Log10 Count')
        clim(ax, [0, global_log_max]); % Apply unified limits (Use 'caxis' if 'clim' fails in 2018b)
    end
end

% --- Figure 3: ZCR / MCR Volatility ---
fig3 = figure('Name', 'Node-Level ZCR Maps', 'Position', [150, 150, 1200, 600], 'Color', 'w');
sgtitle('ZCR/MCR Volatility Features', 'FontSize', 18, 'FontWeight', 'bold');
plot_row_on_surf(1, 2, NodeFeatures.ZCR_Mean, sr, sl, gr, gl, {'Global Destabilizer', 'Score (Mean MCR)'}, parula);
plot_row_on_surf(2, 2, NodeFeatures.ZCR_Top5, sr, sl, gr, gl, {'Targeted Destabilizer', '(Top-5 MCR Mean)'}, parula);

% --- Figure 4: Non-Linear Complexity ---
fig4 = figure('Name', 'Node-Level Polynomial Maps', 'Position', [200, 200, 1200, 600], 'Color', 'w');
sgtitle('Polynomial / Non-Linear Complexity', 'FontSize', 18, 'FontWeight', 'bold');
plot_row_on_surf(1, 2, NodeFeatures.Poly_Footprint, sr, sl, gr, gl, {'Non-Linear Footprint', '(Mean Prevalence)'}, hot);
plot_row_on_surf(2, 2, NodeFeatures.Poly_Entropy, sr, sl, gr, gl, {'Complexity Entropy', '(Diversity of Degrees)'}, turbo);


%% --- HELPER FUNCTIONS ---

function plot_row_on_surf(row_idx, total_rows, plot_vals, sr, sl, gr, gl, tit, cmap)
    % This function places 4 brain views into a specific row of a subplot grid
    
    % Map 200 ROIs to vertices
    cr = zeros(size(gr.cdata)); cl = zeros(size(gl.cdata));
    cr(gr.cdata ~= 0) = plot_vals(gr.cdata(gr.cdata ~= 0));
    cl(gl.cdata ~= 0) = plot_vals(gl.cdata(gl.cdata ~= 0));
    
    % Determine color limits
    c_lims = [min(plot_vals), max(plot_vals)];
    if diff(c_lims) == 0, c_lims = [0 1]; end
    if contains(tit, 'Polarity'), c_lims = [-1 1]; end % Force symmetric scale for polarity
    
    is_gifti = isfield(sr, 'data');
    
    % Calculate Subplot Indices for this row (4 columns per row)
    idx_base = (row_idx - 1) * 4;
    
    % View 1: Left Lateral
    ax1 = subplot(total_rows, 4, idx_base + 1);
    plot_surf(sl, cl, is_gifti, -90, c_lims);
    colormap(ax1, cmap); % FIX: Apply explicitly to axis 1
    
    % Add Row Title elegantly to the left of the first brain
    text(ax1, -0.4, 0.5, tit, ...
        'Units', 'normalized', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
    
    % View 2: Left Medial
    ax2 = subplot(total_rows, 4, idx_base + 2);
    plot_surf(sl, cl, is_gifti, 90, c_lims);
    colormap(ax2, cmap); % FIX: Apply explicitly to axis 2
    
    % View 3: Right Lateral
    ax3 = subplot(total_rows, 4, idx_base + 3);
    plot_surf(sr, cr, is_gifti, 90, c_lims);
    colormap(ax3, cmap); % FIX: Apply explicitly to axis 3
    
    % View 4: Right Medial
    ax4 = subplot(total_rows, 4, idx_base + 4);
    plot_surf(sr, cr, is_gifti, -90, c_lims);
    colormap(ax4, cmap); % FIX: Apply explicitly to axis 4
    
    % Add colorbar next to the last brain in the row
    c = colorbar(ax4, 'Position', [0.92, ax4.Position(2)+0.05, 0.01, ax4.Position(4)-0.1]);
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