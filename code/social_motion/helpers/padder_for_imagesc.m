function aligned_mat = padder_for_imagesc(sorted_cells, sorted_lengths, direction)

% Get the new maximum length (optional, usually same as before)
max_len = max(sorted_lengths);

% Prepare padded matrix
aligned_mat = NaN(length(sorted_lengths), max_len);

% Create indices
row_idx = repelem((1:length(sorted_lengths))', sorted_lengths);

if strcmp(direction, 'offset')
    col_idx = arrayfun(@(len) max_len - len + 1 : max_len, sorted_lengths, 'UniformOutput', false);
elseif strcmp(direction, 'onset')
    col_idx = arrayfun(@(len) 1: len, sorted_lengths, 'UniformOutput', false);
end

col_idx = horzcat(col_idx{:})';

% Concatenate values from sorted cells
all_vals = vertcat(sorted_cells{:});

% Assign values
lin_idx = sub2ind(size(aligned_mat), row_idx, col_idx);
aligned_mat(lin_idx) = all_vals;