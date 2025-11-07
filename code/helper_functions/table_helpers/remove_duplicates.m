function bouts = remove_duplicates(bouts)

[~, idx] = unique(bouts(:, {'fly', 'nloom', 'durations'}), 'rows', 'stable');

% Keep only unique rows
bouts = bouts(idx, :);