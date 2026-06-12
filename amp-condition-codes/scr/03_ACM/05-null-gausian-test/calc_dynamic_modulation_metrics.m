%% Suggested Name: calc_dynamic_modulation_metrics.m
% Calculates Amplitude-based Dynamic Modulation Metrics for fMRI data.
% Supports selecting specific metrics to drastically reduce computation time.
% 
% Usage Examples:
% metrics = calc_dynamic_modulation_metrics(ts, 'n_levels', 20, 'req_metrics', {'variance'});
% metrics = calc_dynamic_modulation_metrics(ts, 'n_levels', 20, 'req_metrics', {'variance', 'zcr'});
% metrics = calc_dynamic_modulation_metrics(ts, 'req_metrics', 'all');

function metrics = calc_dynamic_modulation_metrics(ts, varargin)
% CALC_DYNAMIC_MODULATION_METRICS (Optimized for Speed)

%% 1. Input Parsing
p = inputParser;
addRequired(p, 'ts', @isnumeric);
addParameter(p, 'n_levels', 20, @isscalar);
addParameter(p, 'req_metrics', {'all'}, @(x) iscell(x) || ischar(x));
parse(p, ts, varargin{:});

n_levels = p.Results.n_levels;
req_metrics = p.Results.req_metrics;

if ischar(req_metrics), req_metrics = {req_metrics}; end
if any(strcmpi(req_metrics, 'all'))
    req_metrics = {'variance', 'zcr', 'slope', 'poly'};
end

do_var   = any(strcmpi(req_metrics, 'variance'));
do_zcr   = any(strcmpi(req_metrics, 'zcr'));
do_slope = any(strcmpi(req_metrics, 'slope'));
do_poly  = any(strcmpi(req_metrics, 'poly'));

do_fitting = do_slope || do_poly;

[ntime, nnodes] = size(ts);
nedges = nnodes * (nnodes - 1) / 2;
mask = triu(true(nnodes), 1);

%% 2. Create Mask for Trivial Edges
[u_map, v_map] = find(mask);
edge_mask = false(nedges, nnodes);
for m = 1:nnodes
    edge_mask(:, m) = (u_map == m) | (v_map == m);
end

%% 3. Conditional Pre-computation for Polynomial Fitting
degrees_to_test = [];
V = cell(1, 4);
inv_V = cell(1, 4);

if do_fitting
    x_bins = (1:n_levels)';
    x_norm = (x_bins - mean(x_bins)) / std(x_bins); 
    
    if do_poly
        degrees_to_test = 1:4; % Degrees 0, 1, 2, 3
    elseif do_slope
        degrees_to_test = 2;   % Only degree 1 (linear)
    end
    
    for d_idx = degrees_to_test
        d = d_idx - 1;
        V{d_idx} = x_norm.^(d:-1:0); 
        inv_V{d_idx} = pinv(V{d_idx}'); % Transposed for fast matrix math
    end
end

%% 4. Unconditional Pre-allocation (Safe for parfor)
% Allocating 16MB arrays takes 0.001 seconds. Safe for parfor slicing.
res_var   = zeros(nedges, nnodes, 'single');
res_zcr   = zeros(nedges, nnodes, 'single');
res_slope = zeros(nedges, nnodes, 'single');
res_poly  = zeros(nedges, nnodes, 'uint8'); 

%% 5. Parallel Loop Over Modulators
parfor m_node = 1:nnodes
    % --- A. Quantization & FC Calculation ---
    local_ts = ts; 
    [~, sort_idx] = sort(local_ts(:, m_node));
    bins = round(linspace(1, ntime + 1, n_levels + 1));
    
    level_corrs = zeros(nedges, n_levels, 'single');
    for L = 1:n_levels
        idx_range = sort_idx(bins(L) : bins(L+1)-1);
        if length(idx_range) > 5
             c_mat = corr(local_ts(idx_range, :));
             level_corrs(:, L) = c_mat(mask); 
        else
             level_corrs(:, L) = NaN;
        end
    end
    
    triv_idx = edge_mask(:, m_node); 
    
    % --- B. Metric 1: Variance ---
    if do_var
        var_val = var(level_corrs, 0, 2, 'omitnan');
        var_val(triv_idx) = NaN;
        res_var(:, m_node) = var_val;
    end
    
    % --- C. Metric 2: Zero-Crossing Rate ---
    if do_zcr
        S = sign(level_corrs);
        S(S == 0) = 1; 
        diff_S = diff(S, 1, 2);
        zcr_val = sum(abs(diff_S) == 2, 2); 
        zcr_val(triv_idx) = NaN;
        res_zcr(:, m_node) = zcr_val;
    end
    
    % --- D. Metric 3 & 4: Fitting (Slope / Polynomial) ---
    if do_fitting
        % ONLY run heavy fillmissing if there is actual missing data
        if any(isnan(level_corrs(:)))
            level_corrs = fillmissing(level_corrs, 'linear', 2);
            level_corrs = fillmissing(level_corrs, 'nearest', 2);
        end
        
        bic_mat = zeros(nedges, 4, 'single');
        slope_val = zeros(nedges, 1, 'single');
        poly_val = zeros(nedges, 1, 'uint8');
        
        for d_idx = degrees_to_test
            d = d_idx - 1;
            betas = level_corrs * inv_V{d_idx};
            
            if d == 1
                slope_val = betas(:, 1); 
            end
            
            if do_poly
                Y_hat = betas * V{d_idx}';
                RSS = sum((level_corrs - Y_hat).^2, 2);
                bic_mat(:, d_idx) = n_levels * log((RSS + 1e-10) / n_levels) + (d + 1) * log(n_levels);
            end
        end
        
        if do_slope
            slope_val(triv_idx) = NaN;
            res_slope(:, m_node) = slope_val;
        end
        
        if do_poly
            % Find min BIC considering ONLY the columns we tested
            valid_bics = bic_mat(:, degrees_to_test);
            [~, best_local_idx] = min(valid_bics, [], 2);
            
            % Map back to degree
            actual_d_idx = degrees_to_test(best_local_idx);
            poly_val = uint8(actual_d_idx - 1);
            
            poly_val(triv_idx) = 255;
            res_poly(:, m_node) = poly_val;
        end
    end
end

%% 6. Pack Outputs into Struct
metrics = struct();
if do_var,   metrics.variance = res_var; end
if do_slope, metrics.slope = res_slope; end
if do_zcr,   metrics.zcr = res_zcr; end
if do_poly,  metrics.best_poly_degree = res_poly; end

end