function [bouts_filtered, mask] = quantilizer_v2(bouts, varargin)
% Parse inputs
opt = inputParser;
addParameter(opt, 'mode', 'global');
addParameter(opt, 'total_quantiles', struct('sm', 3, 'fs', 2));
addParameter(opt, 'indexed_quantile', struct('ls', 1, 'fs', 1, 'sm', 1));
addParameter(opt, 'edges', struct()); % Initialize as empty struct
parse(opt, varargin{:});

q = opt.Results.total_quantiles;
target = opt.Results.indexed_quantile;
inputted_edges = opt.Results.edges;

% 1. Create the Reference Table (Only used if edges aren't provided)
% 1. Create the Reference Table
if strcmpi(opt.Results.mode, 'global') || ~isfield(target,'ls')
    ref = bouts;
else
    ref = bouts(bouts.ls == target.ls, :);
end

% 2. Generate or Assign Bin Edges
% Check if fs edges were provided in the struct
if isfield(inputted_edges, 'fs')
    edges_fs = inputted_edges.fs;
else
    edges_fs = quantile(ref.fs, linspace(0, 1, q.fs + 1));
end

% Check if sm edges were provided in the struct
if isfield(inputted_edges, 'sm')
    edges_sm = inputted_edges.sm;
else
    edges_sm = quantile(ref.sm, linspace(0, 1, q.sm + 1));
end

% Handle nloom edges
if isfield(target, 'ln')
    if isfield(inputted_edges, 'ln')
        edges_ln = inputted_edges.ln;
    else
        edges_ln = 0:5:20; % Your default fallback
    end
end

% 3. Tag the Original Table
% Note: values outside edges will result in NaN bins
bouts.bin_fs = discretize(bouts.fs, edges_fs);
bouts.bin_sm = discretize(bouts.sm, edges_sm);

% 4. Master Filter Logic

if isfield(target,'ls')
    is_target_ls = (bouts.ls == target.ls);
else
    is_target_ls = true(height(bouts),1);
end

is_target_fs = (bouts.bin_fs == target.fs);
is_target_sm = (bouts.bin_sm == target.sm);

final_mask = is_target_ls & is_target_fs & is_target_sm;

% Handle nloom if requested
if isfield(target, 'ln')
    bouts.bin_ln = discretize(bouts.nloom, edges_ln);
    final_mask = final_mask & (bouts.bin_ln == target.ln);
end

% 5. Final Output
bouts_filtered = bouts(final_mask, :);
mask = find(final_mask);
end