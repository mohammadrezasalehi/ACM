function summary = get_top_stats(score_mat, top_k)
    % GET_TOP_STATS Extracts summary statistics for top-k candidates.
    % Input: score_mat [edges x nodes], top_k (int)
    % Output: struct with .mean_score, .max_score, .top_indices
    
    [sorted_vals, sorted_idx] = sort(score_mat(:), 'descend');
    
    summary.max_score = sorted_vals(1);
    summary.mean_top = mean(sorted_vals(1:top_k));
    summary.top_indices = sorted_idx(1:top_k); % Store linear indices for Jaccard
end