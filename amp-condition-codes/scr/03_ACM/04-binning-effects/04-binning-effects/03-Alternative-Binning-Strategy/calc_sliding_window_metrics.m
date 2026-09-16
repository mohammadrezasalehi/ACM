%% Suggested Name: calc_sliding_window_metrics.m
% Alternative Binning Strategy: Sliding Window on Amplitude-Sorted Data.
% Generates a high-resolution, smoothed modulation spectrum.

function metrics = calc_sliding_window_metrics(ts, varargin)

%% 1. Input Parsing
p = inputParser;
addRequired(p, 'ts', @isnumeric);
addParameter(p, 'window_size', 120, @isscalar); % Default: 10% of 1200 TRs
addParameter(p, 'step_size', 20, @isscalar);    % Default: slide by 20 TRs
addParameter(p, 'req_metrics', {'all'}, @(x) iscell(x) || ischar(x));
parse(p, ts, varargin{:});

w_size = p.Results.window_size;
s_size = p.Results.step_size;
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

%% 2. Calculate Number of Windows
% Formula for sliding window count: W = floor((N - window_size) / step_size) + 1
n_windows = floor((ntime - w_size) / s_size) + 1;

%% 3. Trivial Edge Mask
[u_map, v_map] = find(mask);
edge_mask = false(nedges, nnodes);
for m = 1:nnodes
    edge_mask(:, m) = (u_map == m) | (v_map == m);
end

%% 4. Polynomial Pre-computation
degrees_to_test = [];
V = cell(1, 4);
inv_V = cell(1, 4);

if do_fitting
    x_bins = (1:n_windows)';
    x_norm = (x_bins - mean(x_bins)) / std(x_bins); 
    
    if do_poly, degrees_to_test = 1:4; 
    elseif do_slope, degrees_to_test = 2; end
    
    for d_idx = degrees_to_test
        d = d_idx - 1;
        V{d_idx} = x_norm.^(d:-1:0); 
        inv_V{d_idx} = pinv(V{d_idx}'); 
    end
end

%% 5. Pre-allocation
if do_var,   res_var   = zeros(nedges, nnodes, 'single'); end
if do_zcr,   res_zcr   = zeros(nedges, nnodes, 'single'); end
if do_slope, res_slope = zeros(nedges, nnodes, 'single'); end
if do_poly,  res_poly  = zeros(nedges, nnodes, 'uint8');  end

%% 6. Parallel Loop
parfor m_node = 1:nnodes
    local_ts = ts; 
    
    % --- DIFFERENCE IS HERE: Sort, then Sliding Window ---
    % 1. Sort the entire time series based on Modulator's amplitude
    [~, sort_idx] = sort(local_ts(:, m_node));
    sorted_ts = local_ts(sort_idx, :);
    
    level_corrs = zeros(nedges, n_windows, 'single');
    
    % 2. Extract correlation for each overlapping window
    for w = 1:n_windows
        start_idx = 1 + (w - 1) * s_size;
        end_idx   = start_idx + w_size - 1;
        
        c_mat = corr(sorted_ts(start_idx:end_idx, :));
        level_corrs(:, w) = c_mat(mask); 
    end
    
    % Mask trivial edges
    triv_idx = edge_mask(:, m_node);
    
    % --- B. Metric Calculations ---
    if do_var
        var_val = var(level_corrs, 0, 2, 'omitnan');
        var_val(triv_idx) = NaN;
        res_var(:, m_node) = var_val;
    end
    
    if do_zcr
        S = sign(level_corrs);
        S(S == 0) = 1; 
        diff_S = diff(S, 1, 2);
        zcr_val = sum(abs(diff_S) == 2, 2); 
        zcr_val(triv_idx) = NaN;
        res_zcr(:, m_node) = zcr_val;
    end
    
    if do_fitting
        bic_mat = zeros(nedges, 4, 'single');
        slope_val = zeros(nedges, 1, 'single');
        
        for d_idx = degrees_to_test
            d = d_idx - 1;
            betas = level_corrs * inv_V{d_idx};
            if d == 1, slope_val = betas(:, 1); end
            
            if do_poly
                Y_hat = betas * V{d_idx}';
                RSS = sum((level_corrs - Y_hat).^2, 2);
                bic_mat(:, d_idx) = n_windows * log((RSS + 1e-10) / n_windows) + (d + 1) * log(n_windows);
            end
        end
        
        if do_slope
            slope_val(triv_idx) = NaN;
            res_slope(:, m_node) = slope_val;
        end
        
        if do_poly
            [~, best_d_idx] = min(bic_mat(:, degrees_to_test), [], 2);
            actual_best_d_idx = degrees_to_test(best_d_idx); 
            best_degree_val = uint8(actual_best_d_idx - 1);
            best_degree_val(triv_idx) = 255;
            res_poly(:, m_node) = best_degree_val;
        end
    end
end

%% 7. Pack Outputs
metrics = struct();
if do_var,   metrics.variance = res_var; end
if do_slope, metrics.slope = res_slope; end
if do_zcr,   metrics.zcr = res_zcr; end
if do_poly,  metrics.best_poly_degree = res_poly; end

end