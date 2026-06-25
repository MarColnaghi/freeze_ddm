function al = align_peaks(clu, sm_freeze_full, bouts_proc, inizio)
% ALIGN_PEAKS  Peak/linear-fit alignment of the clustered pre-freeze SM, and
% the magnitude-sorted grouping used by the peak-vs-magnitude figures.
%
%   al = align_peaks(clu, sm_freeze_full, bouts_proc, inizio)
%
% clu is the struct returned by cluster_sm (uses sort_order, valid_mask,
% sorted_matrix). Fields of AL (the outputs consumed downstream):
%   valid_sorted, x_axis_peak, aligned_peak_valid, aligned_peak_mag,
%   sorted_matrix_valid, break_times_valid, break_times_peak_valid,
%   break_times_peak_mag, n_rows_valid, n_rows_mag, n_cols, n_mag_clusters,
%   idx_mag_cluster_sorted, boundaries_sm

sort_order    = clu.sort_order;
valid_mask    = clu.valid_mask;
sorted_matrix = clu.sorted_matrix;

% Individual peaks (shifted to the plot's x-start)
[~, peak_idx_raw] = max(sm_freeze_full(sort_order, :), [], 2);
individual_peaks = inizio + (peak_idx_raw - 1);
peak_buffer = 30; %#ok<NASGU>

valid_sorted = valid_mask(sort_order);   % already in sorted order

% --- 2. -------- PEAK ALIGNMENT --------
ref_peak = 315;
shifts_peak = round(ref_peak - individual_peaks);

% --- 3. -------- LINEAR FIT ALIGNMENT --------
x = individual_peaks(valid_sorted);
y = find(valid_sorted);

coeffs = polyfit(x, y, 1);
slope = coeffs(1);
intercept = coeffs(2);

n_rows = size(sorted_matrix,1);
n_cols = size(sorted_matrix,2);

y_all = (1:n_rows)';
x_fit_all = (y_all - intercept) / slope;

ref_fit = mean(x_fit_all(valid_sorted));
shifts_fit = round(ref_fit - x_fit_all);

% --- 4. -------- PAD MATRICES (NO CROPPING) --------
max_shift = ceil(max(abs([shifts_peak(:); shifts_fit(:)]), [], 'omitnan'));
pad = max_shift + 30;

new_cols = n_cols + 2*pad;

aligned_peak = nan(n_rows, new_cols);
aligned_fit  = nan(n_rows, new_cols);

base_idx = pad + (1:n_cols);

for i = 1:n_rows
    if valid_sorted(i)
        if ~isnan(shifts_peak(i))
            aligned_peak(i, base_idx + shifts_peak(i)) = sorted_matrix(i,:);
        end

        if ~isnan(shifts_fit(i))
            aligned_fit(i, base_idx + shifts_fit(i)) = sorted_matrix(i,:);
        end
    end
end

% --- 5. -------- CENTERED TIME AXIS --------
ref_col_peak = pad + round(ref_peak);
ref_col_fit  = pad + round(ref_fit); %#ok<NASGU>

x_axis_peak = (1:new_cols) - ref_col_peak;

% --- 6. -------- BREAK TIMES (CENTERED) --------
break_times = bouts_proc.durations(sort_order);
break_times(break_times > 630) = NaN;

break_times_peak = break_times + shifts_peak - round(ref_peak);

% --- 7. -------- FILTER VALID TRIALS --------
sorted_matrix_valid = sorted_matrix(valid_sorted, :);
aligned_peak_valid  = aligned_peak(valid_sorted, :);

break_times_valid       = break_times(valid_sorted);
break_times_peak_valid  = break_times_peak(valid_sorted);

n_rows_valid = sum(valid_sorted);

% --- 10. -------- SORT BY PEAK MAGNITUDE --------
peak_mag = max(aligned_peak_valid, [], 2, 'omitnan');
[peak_mag_sorted, sort_idx_mag] = sort(peak_mag, 'descend');

aligned_peak_mag = aligned_peak_valid(sort_idx_mag, :);
break_times_peak_mag = break_times_peak_valid(sort_idx_mag);

n_rows_mag = size(aligned_peak_mag,1);

n_mag_clusters = 6;
N = numel(peak_mag_sorted);
group_size = ceil(N / n_mag_clusters);

idx_mag_cluster_sorted = zeros(N,1);

for k = 1:n_mag_clusters
    inds = (k-1)*group_size + 1 : min(k*group_size, N);
    idx_mag_cluster_sorted(inds) = k;
end

edges = find(diff(idx_mag_cluster_sorted) ~= 0);
boundaries_sm = [0; edges; N];

al = struct('valid_sorted', valid_sorted, 'x_axis_peak', x_axis_peak, ...
    'aligned_peak_valid', aligned_peak_valid, 'aligned_peak_mag', aligned_peak_mag, ...
    'sorted_matrix_valid', sorted_matrix_valid, ...
    'break_times_valid', break_times_valid, 'break_times_peak_valid', break_times_peak_valid, ...
    'break_times_peak_mag', break_times_peak_mag, ...
    'n_rows_valid', n_rows_valid, 'n_rows_mag', n_rows_mag, 'n_cols', n_cols, ...
    'n_mag_clusters', n_mag_clusters, 'idx_mag_cluster_sorted', idx_mag_cluster_sorted, ...
    'boundaries_sm', boundaries_sm);
end
