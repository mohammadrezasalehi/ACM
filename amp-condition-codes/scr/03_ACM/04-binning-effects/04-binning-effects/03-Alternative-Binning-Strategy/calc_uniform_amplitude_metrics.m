%% Suggested Name: calc_uniform_amplitude_metrics.m
% Alternative Binning Strategy: Uniform Amplitude (Value-based) Binning.
% Features custom ultra-fast interpolation for empty bins.

function metrics = calc_uniform_amplitude_metrics(ts, varargin)

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

%% 2. Trivial Edge Mask
[u_map, v_map] = find(mask);
edge_mask = false(nedges, nnodes);
for m = 1:nnodes
    edge_mask(:, m) = (u_map == m) | (v_map == m);
end

%% 3. Polynomial Pre-computation
degrees_to_test = [];
V = cell(1, 4);
inv_V = cell(1, 4);

if do_fitting
    x_bins = (1:n_levels)';
    x_norm = (x_bins - mean(x_bins)) / std(x_bins); 
    if do_poly, degrees_to_test = 1:4; 
    elseif do_slope, degrees_to_test = 2; end
    
    for d_idx = degrees_to_test
        d = d_idx - 1;
        V{d_idx} = x_norm.^(d:-1:0); 
        inv_V{d_idx} = pinv(V{d_idx}'); 
    end
end

%% 4. Pre-allocation
if do_var,   res_var   = zeros(nedges, nnodes, 'single'); end
if do_zcr,   res_zcr   = zeros(nedges, nnodes, 'single'); end
if do_slope, res_slope = zeros(nedges, nnodes, 'single'); end
if do_poly,  res_poly  = zeros(nedges, nnodes, 'uint8');  end

% Array to keep track of NaN percentage for each modulator
nan_percentages = zeros(nnodes, 1);

%% 5. Parallel Loop
parfor m_node = 1:nnodes
    local_ts = ts; 
    mod_signal = local_ts(:, m_node);
    
    % --- Uniform Amplitude Binning ---
    min_val = min(mod_signal);
    max_val = max(mod_signal);
    edges = linspace(min_val, max_val, n_levels + 1);
    edges(1) = -Inf; edges(end) = Inf;
    bin_assignments = discretize(mod_signal, edges);
    
    level_corrs = zeros(nedges, n_levels, 'single');
    for L = 1:n_levels
        idx_range = find(bin_assignments == L);
        
        if length(idx_range) > 2
             c_mat = corr(local_ts(idx_range, :));
             level_corrs(:, L) = c_mat(mask); 
        else
             level_corrs(:, L) = NaN; % Empty or sparse bin
        end
    end
    
    % --- Calculate NaN Percentage before interpolation ---
    total_elements = nedges * n_levels;
    nan_count = sum(isnan(level_corrs(:)));
    nan_percentages(m_node) = (nan_count / total_elements) * 100;
    
    % --- ULTRA-FAST CUSTOM INTERPOLATION (Vectorized) ---
    % Find missing values
    missing_mask = isnan(level_corrs);
    
    % Iterate only a few times to fill gaps (usually 1 or 2 iterations is enough)
    % This logic averages the left and right neighbors for NaNs
    for iter = 1:3 
        if ~any(missing_mask(:)), break; end
        
        % Shift matrices left and right (padding with NaN at edges)
        left_neighbor = [NaN(nedges, 1), level_corrs(:, 1:end-1)];
        right_neighbor = [level_corrs(:, 2:end), NaN(nedges, 1)];
        
        % Compute mean of available neighbors (omitnan automatically ignores missing ones)
        % We stack them into a 3D array [nedges x n_levels x 2] and take mean along 3rd dim
        neighbors = cat(3, left_neighbor, right_neighbor);
        mean_neighbors = mean(neighbors, 3, 'omitnan');
        
        % Only update the originally missing values that now have a valid neighbor mean
        update_mask = missing_mask & ~isnan(mean_neighbors);
        level_corrs(update_mask) = mean_neighbors(update_mask);
        
        % Update missing mask for next iteration
        missing_mask = isnan(level_corrs);
    end
    
    % If any NaNs remain (e.g. an entire row was empty), replace with 0 to prevent downstream failure
    level_corrs(isnan(level_corrs)) = 0; 
    
    % --- Apply Trivial Edge Mask ---
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
                bic_mat(:, d_idx) = n_levels * log((RSS + 1e-10) / n_levels) + (d + 1) * log(n_levels);
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

%% 6. Pack Outputs
metrics = struct();
if do_var,   metrics.variance = res_var; end
if do_slope, metrics.slope = res_slope; end
if do_zcr,   metrics.zcr = res_zcr; end
if do_poly,  metrics.best_poly_degree = res_poly; end

% Attach NaN diagnostic info to the output struct
metrics.nan_percent_per_node = nan_percentages;

end