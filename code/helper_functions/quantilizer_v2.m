function [bouts_filtered, mask] = quantilizer_v2(bouts, varargin)
    
% Parse inputs
    opt = inputParser;
    addParameter(opt, 'mode', 'global'); 
    addParameter(opt, 'total_quantiles', struct('sm', 3, 'fs', 2)); 
    addParameter(opt, 'indexed_quantile', struct('ls', 1, 'fs', 1, 'sm', 1));
    parse(opt, varargin{:});
    
    q = opt.Results.total_quantiles;
    target = opt.Results.indexed_quantile;

    % 1. Create the Reference Table
    if strcmpi(opt.Results.mode, 'global')
        ref = bouts;
    else
        ref = bouts(bouts.ls == target.ls, :);
    end

    % 2. Generate Bin Edges (The "Slicers")
    edges_fs = quantile(ref.fs, linspace(0, 1, q.fs + 1));
    edges_sm = quantile(ref.sm, linspace(0, 1, q.sm + 1));
    
    % For nloom, we'll use your custom 0:5:20 or dynamic quantiles
    if isfield(target, 'ln')
        edges_ln = 0:5:20; 
    end

    % 3. Tag the Original Table
    % This adds new columns to your table identifying which bin each row is in
    bouts.bin_fs = discretize(bouts.fs, edges_fs);
    bouts.bin_sm = discretize(bouts.sm, edges_sm);
    
    % 4. Master Filter Logic
    % Now filtering is just a simple table query
    is_target_loom = (bouts.ls == target.ls);
    is_target_fs   = (bouts.bin_fs == target.fs);
    is_target_sm   = (bouts.bin_sm == target.sm);
    
    final_mask = is_target_loom & is_target_fs & is_target_sm;
    
    % Handle nloom if requested
    if isfield(target, 'ln')
        bouts.bin_ln = discretize(bouts.nloom, edges_ln);
        final_mask = final_mask & (bouts.bin_ln == target.ln);
    end

    % 5. Final Output
    bouts_filtered = bouts(final_mask, :);
    mask = find(final_mask);
end