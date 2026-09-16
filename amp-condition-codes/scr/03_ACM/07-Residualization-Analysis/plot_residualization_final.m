%% =========================================================================
% plot_residualization_final.m
%
% PURPOSE:
%   Final publication-ready visualization of the residualization analysis.
%
%   PANEL B:
%       Variance explained by:
%           1) Static FC predictors
%           2) Static FC + AR(1) predictors
%           3) Incremental Delta R^2 attributable to AR(1)
%
%   PANEL C:
%       Preservation of the original spatial topology after residualization:
%           Pearson r between Original and Residualized maps at:
%           1) Raw triplet level
%           2) Node level
%           3) Network level
%
% IMPORTANT:
%   This script ONLY READS SAVED RESULTS.
%   It DOES NOT rerun any regression/residualization.
%
% MATLAB:
%   Compatible with MATLAB R2018b
%
% =========================================================================

clc;
clear;
close all;

fprintf('\n============================================================\n');
fprintf(' FINAL RESIDUALIZATION VISUALIZATION\n');
fprintf('============================================================\n\n');


%% =========================================================================
% 1. CONFIGURATION
% =========================================================================

project_root = 'F:\PhD Code\My_PhD_Project';

session   = 'REST1_LR';
clean_pipe = 'No_clean\';

n_rois  = 200;
n_nets  = 17;

metrics = {'variance', 'slope', 'zcr', 'best_poly_degree'};

metric_labels = { ...
    'Variance', ...
    'Slope', ...
    'Volatility (ZCR)', ...
    'Nonlinear Complexity'};

%% -------------------------------------------------------------------------
% Directories
% -------------------------------------------------------------------------

orig_group_dir = fullfile( ...
    project_root, ...
    'results', ...
    'Modulated_FC', ...
    'Group_Level', ...
    session);

resid_dir = fullfile( ...
    project_root, ...
    'results', ...
    'Modulated_FC_Residuals', ...
    session, ...
    clean_pipe);

r2_file = fullfile( ...
    project_root, ...
    'results', ...
    'Modulated_FC_Residuals', ...
    'R2_Summary_Statistics.mat');

%% -------------------------------------------------------------------------
% Atlas information
% -------------------------------------------------------------------------

atlas_info_path = fullfile( ...
    project_root, ...
    'data', ...
    'atlases', ...
    'Schaefer2018_200Parcels_Kong2022_17Networks_order_info.txt');

utility_path = fullfile( ...
    project_root, ...
    'src', ...
    '04_utility');

if exist(utility_path, 'dir')
    addpath(genpath(utility_path));
end

%% -------------------------------------------------------------------------
% Output directory
% -------------------------------------------------------------------------

output_dir = fullfile( ...
    project_root, ...
    'results', ...
    'Modulated_FC_Residuals', ...
    'Final_Figure');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


%% =========================================================================
% 2. BASIC CHECKS
% =========================================================================

fprintf('Checking required files...\n');

if ~exist(r2_file, 'file')
    error('R2 summary file not found:\n%s', r2_file);
end

if ~exist(orig_group_dir, 'dir')
    error('Original group directory not found:\n%s', orig_group_dir);
end

if ~exist(resid_dir, 'dir')
    error('Residual directory not found:\n%s', resid_dir);
end

fprintf('  R2 file       : OK\n');
fprintf('  Original data : OK\n');
fprintf('  Residual data : OK\n\n');


%% =========================================================================
% 3. LOAD R2 SUMMARY
% =========================================================================

fprintf('Loading R2 summary...\n');

R2 = load(r2_file);

%% -------------------------------------------------------------------------
% Calculate mean R2 values across subjects
% -------------------------------------------------------------------------

mean_FC    = nan(1,4);
mean_Full  = nan(1,4);
mean_Delta = nan(1,4);

for m = 1:4

    metric_name = metrics{m};

    if ~isfield(R2, metric_name)
        error('Metric "%s" was not found in R2 summary.', metric_name);
    end

    S = R2.(metric_name);

    %% Static FC R2
    if isfield(S, 'R2_FC')
        mean_FC(m) = mean(S.R2_FC(~isnan(S.R2_FC)));
    else
        error('R2_FC not found for metric %s.', metric_name);
    end

    %% Full model
    if isfield(S, 'R2_Full')

        mean_Full(m) = mean(S.R2_Full(~isnan(S.R2_Full)));

    elseif isfield(S, 'R2_Full_Model')

        mean_Full(m) = mean(S.R2_Full_Model(~isnan(S.R2_Full_Model)));

    elseif isfield(S, 'Delta_R2')

        mean_Delta(m) = mean(S.Delta_R2(~isnan(S.Delta_R2)));
        mean_Full(m) = mean_FC(m) + mean_Delta(m);

    elseif isfield(S, 'R2_AR')

        % In the existing residualization pipeline R2_AR represents
        % the incremental contribution of AR(1).
        mean_Delta(m) = mean(S.R2_AR(~isnan(S.R2_AR)));
        mean_Full(m) = mean_FC(m) + mean_Delta(m);

    else
        error(['Could not identify the full-model R2 or Delta R2 for ', ...
               'metric %s.'], metric_name);
    end

    %% If Delta was not explicitly calculated yet
    if isnan(mean_Delta(m))
        mean_Delta(m) = mean_Full(m) - mean_FC(m);
    end

end

%% Convert to percentage

Data_Panel_B = [ ...
    mean_FC(:), ...
    mean_Full(:), ...
    mean_Delta(:)] * 100;

fprintf('\nVariance explained:\n');

for m = 1:4
    fprintf('%-22s FC = %6.2f%% | Full = %6.2f%% | Delta = %6.2f%%\n', ...
        metric_labels{m}, ...
        Data_Panel_B(m,1), ...
        Data_Panel_B(m,2), ...
        Data_Panel_B(m,3));
end


%% =========================================================================
% 4. LOAD ORIGINAL GROUP-LEVEL MAPS
% =========================================================================

fprintf('\nLoading original group-level maps...\n');

Orig_V = load(fullfile( ...
    orig_group_dir, ...
    'Group_Variance_Consensus.mat'));

Orig_S = load(fullfile( ...
    orig_group_dir, ...
    'Group_Slope_Tmap.mat'));

Orig_Z = load(fullfile( ...
    orig_group_dir, ...
    'Group_ZCR_Median.mat'));

Orig_P = load(fullfile( ...
    orig_group_dir, ...
    'Group_Poly_Prevalence.mat'));


%% -------------------------------------------------------------------------
% Extract matrices
%
% Each matrix should be:
%       19900 x 200
%
% rows    = target edges
% columns = modulators
% -------------------------------------------------------------------------

Orig_M = cell(1,4);

Orig_M{1} = double(Orig_V.Consensus_Percent);
Orig_M{2} = double(Orig_S.T_map);
Orig_M{3} = double(Orig_Z.Median);

%% For polynomial map, use Mode as in the existing group-level pipeline
if isfield(Orig_P, 'Mode')
    Orig_M{4} = double(Orig_P.Mode);
elseif isfield(Orig_P, 'Prevalence_Percent')
    Orig_M{4} = double(Orig_P.Prevalence_Percent);
else
    error('Polynomial group-level variable not found.');
end


%% =========================================================================
% 5. LOAD AND AGGREGATE RESIDUALIZED SUBJECT MAPS
% =========================================================================

fprintf('\nLoading residualized subject maps...\n');

resid_files = dir(fullfile(resid_dir, '*_residualized.mat'));

n_subs = length(resid_files);

if n_subs == 0
    error('No residualized subject files were found.');
end

fprintf('  Number of subjects: %d\n', n_subs);

%% -------------------------------------------------------------------------
% Number of edges
% -------------------------------------------------------------------------

[u_map, v_map] = find(triu(true(n_rois),1));

n_edges = length(u_map);

fprintf('  Number of ROIs    : %d\n', n_rois);
fprintf('  Number of edges   : %d\n', n_edges);


%% -------------------------------------------------------------------------
% Network assignment
% -------------------------------------------------------------------------

if ~exist(atlas_info_path, 'file')
    error('Atlas information file not found:\n%s', atlas_info_path);
end

[~, roi_nets] = load_atlas_info(atlas_info_path);

unique_nets = unique(roi_nets, 'stable');

if length(unique_nets) ~= n_nets
    warning('Expected %d networks, found %d.', ...
        n_nets, length(unique_nets));
end

net_idx_map = zeros(n_rois,1);

for i = 1:n_rois
    idx = find(strcmp(unique_nets, roi_nets{i}), 1);

    if isempty(idx)
        error('Network assignment missing for ROI %d.', i);
    end

    net_idx_map(i) = idx;
end


%% =========================================================================
% 6. GROUP MEAN OF RESIDUALIZED MATRICES
%
% IMPORTANT:
%   NaNs are NOT replaced by zero.
%   Each matrix element is averaged over the subjects for whom it is valid.
% =========================================================================

Resid_M = cell(1,4);

for m = 1:4

    fprintf('\nAggregating residualized metric: %s\n', ...
        metric_labels{m});

    sum_mat   = zeros(n_edges, n_rois, 'double');
    count_mat = zeros(n_edges, n_rois, 'double');

    for s = 1:n_subs

        R = load(fullfile( ...
            resid_files(s).folder, ...
            resid_files(s).name));

        metric_name = metrics{m};

        if ~isfield(R, metric_name)
            error('Variable "%s" missing in %s.', ...
                metric_name, resid_files(s).name);
        end

        X = double(R.(metric_name));

        if ~isequal(size(X), [n_edges, n_rois])
            error(['Unexpected matrix size for %s in %s. ', ...
                   'Expected %d x %d.'], ...
                metric_name, ...
                resid_files(s).name, ...
                n_edges, ...
                n_rois);
        end

        %% Invalid polynomial code
        if m == 4
            valid = isfinite(X) & X ~= 255;
        else
            valid = isfinite(X);
        end

        sum_mat(valid)   = sum_mat(valid) + X(valid);
        count_mat(valid) = count_mat(valid) + 1;

        clear X R;
    end

    %% NaN-aware group mean
    G = nan(n_edges, n_rois);

    valid_count = count_mat > 0;

    G(valid_count) = ...
        sum_mat(valid_count) ./ count_mat(valid_count);

    Resid_M{m} = G;

    clear sum_mat count_mat G;

    fprintf('  Done.\n');
end


%% =========================================================================
% 7. CALCULATE RESIDUAL PRESERVATION AT THREE LEVELS
%
% LEVEL 1: RAW
%       Full 19900 x 200 triplet-level map
%
% LEVEL 2: NODE
%       200-element modulator profile
%
% LEVEL 3: NETWORK
%       17 x 17 network matrix
% =========================================================================

fprintf('\n============================================================\n');
fprintf(' Calculating Original vs Residualized Correspondence\n');
fprintf('============================================================\n');

Data_Panel_C = nan(4,3);

% Additional storage for diagnostics
Pearson_Raw     = nan(4,1);
Pearson_Node    = nan(4,1);
Pearson_Network = nan(4,1);

for m = 1:4

    fprintf('\nMetric: %s\n', metric_labels{m});

    Orig = Orig_M{m};
    Resid = Resid_M{m};

    %% ---------------------------------------------------------------------
    % INVALID VALUES
    % ---------------------------------------------------------------------

    if m == 4
        valid_raw = ...
            isfinite(Orig) & ...
            isfinite(Resid) & ...
            Orig ~= 255 & ...
            Resid ~= 255;
    else
        valid_raw = ...
            isfinite(Orig) & ...
            isfinite(Resid);
    end

    %% ---------------------------------------------------------------------
    % 1. RAW LEVEL
    %
    % Full 19900 x 200 map
    % ---------------------------------------------------------------------

    x = Orig(valid_raw);
    y = Resid(valid_raw);

    if length(x) < 3
        error('Insufficient valid raw-level observations.');
    end

    r_raw = corr(x, y, 'Type', 'Pearson');

    Pearson_Raw(m) = r_raw;

    fprintf('  Raw      : r = %.4f\n', r_raw);

    clear x y;


    %% ---------------------------------------------------------------------
    % 2. NODE LEVEL
    %
    % Exactly follows the aggregation logic used in your pipeline:
    %
    % Variance          -> SD across edges
    % Slope             -> mean absolute slope
    % ZCR               -> mean
    % Polynomial degree -> mean
    % ---------------------------------------------------------------------

    switch m

        case 1
            Orig_Node = std(Orig, 0, 1, 'omitnan')';
            Resid_Node = std(Resid, 0, 1, 'omitnan')';

        case 2
            Orig_Node = mean(abs(Orig), 1, 'omitnan')';
            Resid_Node = mean(abs(Resid), 1, 'omitnan')';

        case 3
            Orig_Node = mean(Orig, 1, 'omitnan')';
            Resid_Node = mean(Resid, 1, 'omitnan')';

        case 4
            Orig_Node = mean(Orig, 1, 'omitnan')';
            Resid_Node = mean(Resid, 1, 'omitnan')';

    end

    valid_node = ...
        isfinite(Orig_Node) & ...
        isfinite(Resid_Node);

    r_node = corr( ...
        Orig_Node(valid_node), ...
        Resid_Node(valid_node), ...
        'Type', 'Pearson');

    Pearson_Node(m) = r_node;

    fprintf('  Node     : r = %.4f\n', r_node);


    %% ---------------------------------------------------------------------
    % 3. NETWORK LEVEL
    %
    % First average across modulators for each target edge.
    % Then map target edges to the 17 x 17 network matrix.
    % ---------------------------------------------------------------------

    Orig_Edge = mean(Orig, 2, 'omitnan');
    Resid_Edge = mean(Resid, 2, 'omitnan');

    if m == 4
        Orig_Edge(Orig_Edge == 255) = NaN;
        Resid_Edge(Resid_Edge == 255) = NaN;
    end

    Orig_Net = map_edges_to_network( ...
        Orig_Edge, ...
        u_map, ...
        v_map, ...
        net_idx_map, ...
        n_nets);

    Resid_Net = map_edges_to_network( ...
        Resid_Edge, ...
        u_map, ...
        v_map, ...
        net_idx_map, ...
        n_nets);

    valid_net = ...
        isfinite(Orig_Net) & ...
        isfinite(Resid_Net);

    r_net = corr( ...
        Orig_Net(valid_net), ...
        Resid_Net(valid_net), ...
        'Type', 'Pearson');

    Pearson_Network(m) = r_net;

    fprintf('  Network  : r = %.4f\n', r_net);

    %% Store
    Data_Panel_C(m,:) = ...
        [r_raw, r_node, r_net];

    clear Orig Resid Orig_Node Resid_Node
    clear Orig_Edge Resid_Edge Orig_Net Resid_Net
end


%% =========================================================================
% 8. DISPLAY FINAL NUMERICAL TABLES
% =========================================================================

fprintf('\n============================================================\n');
fprintf(' FINAL RESULTS\n');
fprintf('============================================================\n');

fprintf('\nPANEL B: VARIANCE EXPLAINED\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-24s %10s %10s %10s\n', ...
    'Metric', 'FC Only', 'FC+AR1', 'Delta R2');

for m = 1:4
    fprintf('%-24s %9.2f%% %9.2f%% %9.2f%%\n', ...
        metric_labels{m}, ...
        Data_Panel_B(m,1), ...
        Data_Panel_B(m,2), ...
        Data_Panel_B(m,3));
end


fprintf('\nPANEL C: RESIDUAL PRESERVATION\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-24s %10s %10s %10s\n', ...
    'Metric', 'Raw', 'Node', 'Network');

for m = 1:4
    fprintf('%-24s %10.4f %10.4f %10.4f\n', ...
        metric_labels{m}, ...
        Data_Panel_C(m,1), ...
        Data_Panel_C(m,2), ...
        Data_Panel_C(m,3));
end


%% =========================================================================
% 9. SAVE NUMERICAL RESULTS
% =========================================================================

Results = struct();

Results.Panel_B_Variance_Explained = Data_Panel_B;
Results.Panel_C_Residual_Preservation = Data_Panel_C;

Results.metric_labels = metric_labels;

Results.R2_FC_mean = mean_FC;
Results.R2_Full_mean = mean_Full;
Results.Delta_R2_mean = mean_Delta;

Results.Pearson_Raw = Pearson_Raw;
Results.Pearson_Node = Pearson_Node;
Results.Pearson_Network = Pearson_Network;

Results.n_subjects = n_subs;
Results.n_rois = n_rois;
Results.n_edges = n_edges;
Results.n_networks = n_nets;

save(fullfile(output_dir, ...
    'Residualization_Final_Statistics.mat'), ...
    'Results', '-v7.3');

fprintf('\nSaved MAT results.\n');


%% =========================================================================
% 10. SAVE CSV TABLES
% =========================================================================

%% Panel B table

PanelB_Table = [ ...
    Data_Panel_B(:,1), ...
    Data_Panel_B(:,2), ...
    Data_Panel_B(:,3)];

csvwrite(fullfile(output_dir, ...
    'Panel_B_Variance_Explained.csv'), ...
    PanelB_Table);


%% Panel C table

PanelC_Table = [ ...
    Data_Panel_C(:,1), ...
    Data_Panel_C(:,2), ...
    Data_Panel_C(:,3)];

csvwrite(fullfile(output_dir, ...
    'Panel_C_Residual_Preservation.csv'), ...
    PanelC_Table);

fprintf('Saved CSV tables.\n');


%% =========================================================================
% 11. PUBLICATION-QUALITY FIGURE
% =========================================================================

fprintf('\nGenerating final figure...\n');

fig = figure( ...
    'Color', 'w', ...
    'Position', [100 100 1250 560], ...
    'Name', 'Residualization Analysis - Final Figure');

%% -------------------------------------------------------------------------
% PANEL B
% -------------------------------------------------------------------------

ax1 = subplot(1,2,1);

imagesc(Data_Panel_B);

%% Custom sequential colormap
ncolors = 256;

cmap_R2 = [ ...
    ones(ncolors,1), ...
    linspace(1,0.20,ncolors)', ...
    linspace(1,0.20,ncolors)'];

colormap(ax1, cmap_R2);

caxis([0 max(100, max(Data_Panel_B(:))*1.05)]);

set(ax1, ...
    'YDir', 'reverse', ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'FontName', 'Arial', ...
    'FontSize', 10, ...
    'LineWidth', 1);

xticks(1:3);
xticklabels({ ...
    'FC only', ...
    'FC + AR1', ...
    '\DeltaR^2 (AR1)'});

yticks(1:4);
yticklabels(metric_labels);

xlabel('Predictor model', ...
    'FontSize', 11, ...
    'FontWeight', 'normal');

ylabel('Modulation metric', ...
    'FontSize', 11, ...
    'FontWeight', 'normal');

title('B  Variance explained', ...
    'FontSize', 13, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');

axis square;

%% Text labels
for r = 1:4
    for c = 1:3

        val = Data_Panel_B(r,c);

        if val > 0.60 * max(Data_Panel_B(:))
            txt_color = 'w';
        else
            txt_color = 'k';
        end

        text(c, r, ...
            sprintf('%.1f%%', val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'Color', txt_color);

    end
end

cb1 = colorbar('southoutside');

set(cb1, ...
    'FontName', 'Arial', ...
    'FontSize', 9);

ylabel(cb1, 'Variance explained (%)');


%% -------------------------------------------------------------------------
% PANEL C
% -------------------------------------------------------------------------

ax2 = subplot(1,2,2);

imagesc(Data_Panel_C);

%% Blue sequential colormap
cmap_R = [ ...
    linspace(1,0.10,ncolors)', ...
    linspace(1,0.40,ncolors)', ...
    ones(ncolors,1)];

colormap(ax2, cmap_R);

caxis([0 1]);

set(ax2, ...
    'YDir', 'reverse', ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'FontName', 'Arial', ...
    'FontSize', 10, ...
    'LineWidth', 1);

xticks(1:3);
xticklabels({ ...
    'Raw', ...
    'Node', ...
    'Network'});

yticks(1:4);
yticklabels(metric_labels);

xlabel('Spatial representation', ...
    'FontSize', 11);

ylabel('Modulation metric', ...
    'FontSize', 11);

title('C  Residual preservation', ...
    'FontSize', 13, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');

axis square;

%% Text labels

for r = 1:4
    for c = 1:3

        val = Data_Panel_C(r,c);

        if val > 0.60
            txt_color = 'w';
        else
            txt_color = 'k';
        end

        text(c, r, ...
            sprintf('%.2f', val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'Color', txt_color);

    end
end

cb2 = colorbar('southoutside');

set(cb2, ...
    'FontName', 'Arial', ...
    'FontSize', 9);

ylabel(cb2, 'Pearson correlation (r)');


%% =========================================================================
% 12. FINAL FIGURE FORMATTING
% =========================================================================

%% Manually tighten subplot positions
set(ax1, 'Position', [0.08 0.20 0.36 0.67]);
set(ax2, 'Position', [0.56 0.20 0.36 0.67]);

%% Panel labels
annotation(fig, ...
    'textbox', ...
    [0.025 0.91 0.04 0.05], ...
    'String', 'B', ...
    'EdgeColor', 'none', ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

annotation(fig, ...
    'textbox', ...
    [0.505 0.91 0.04 0.05], ...
    'String', 'C', ...
    'EdgeColor', 'none', ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');


%% =========================================================================
% 13. SAVE FIGURE
% =========================================================================

fig_base = fullfile(output_dir, ...
    'Figure3_Residualization_Heatmaps');

%% MATLAB figure
savefig(fig, [fig_base '.fig']);

%% TIFF - publication quality
print(fig, ...
    [fig_base '.tif'], ...
    '-dtiff', ...
    '-r600');

%% PDF / vector
print(fig, ...
    [fig_base '.pdf'], ...
    '-dpdf', ...
    '-painters');

fprintf('\n============================================================\n');
fprintf(' FINAL FIGURE SAVED\n');
fprintf('============================================================\n');
fprintf('%s\n', output_dir);
fprintf('\n');


%% =========================================================================
% LOCAL FUNCTION
% =========================================================================

function net_mat = map_edges_to_network( ...
    edge_vec, ...
    u, ...
    v, ...
    net_idx_map, ...
    n_nets)

% -------------------------------------------------------------------------
% Maps a 19900 x 1 edge vector into a 17 x 17 network matrix.
%
% The implementation follows the network aggregation logic used in the
% existing robustness pipeline:
%
%   - each edge contributes to its corresponding network pair
%   - between-network values are mirrored
%   - values are averaged within each network pair
% -------------------------------------------------------------------------

sum_scores   = zeros(n_nets, n_nets);
count_scores = zeros(n_nets, n_nets);

for e = 1:length(u)

    ni = u(e);
    nj = v(e);

    net_i = net_idx_map(ni);
    net_j = net_idx_map(nj);

    val = edge_vec(e);

    if ~isfinite(val)
        continue;
    end

    %% i-j
    sum_scores(net_i, net_j) = ...
        sum_scores(net_i, net_j) + val;

    count_scores(net_i, net_j) = ...
        count_scores(net_i, net_j) + 1;

    %% j-i
    if net_i ~= net_j

        sum_scores(net_j, net_i) = ...
            sum_scores(net_j, net_i) + val;

        count_scores(net_j, net_i) = ...
            count_scores(net_j, net_i) + 1;

    end

end

net_mat = nan(n_nets, n_nets);

valid = count_scores > 0;

net_mat(valid) = ...
    sum_scores(valid) ./ count_scores(valid);

end