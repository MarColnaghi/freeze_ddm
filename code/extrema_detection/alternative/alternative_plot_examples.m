
% Load the table with the computed log likelihoods
model_list = {'dddim2', 'ded2'};
run_list   = {'run02_260325', 'run05_260325'};
version = {''};

paths_analysis = path_generator('folder', fullfile('extrema_detection','alternative_test', model_list{2}, strcat(run_list{2}, version{1})));
paths_analysis.fig = fullfile(paths_analysis.fig, 'examples');
col = cmapper('', 2);

ll = struct2table(importdata(fullfile(paths_analysis.results, 'll_table.mat')));
data = alternative_test_collect_data(model_list, run_list);

% --- 3. Main Loop: Iterate over Bouts ---
freezes_ref = data(1).table;
n_bouts     = height(freezes_ref);
censoring_val = data(1).results.points.censoring;
res_points = data(1).results.points; % Slice to avoid passing large 'data' struct

% Slice the tables to reduce communication overhead
out_all_tv = data(1).out_all_tv;
out_all_st = data(1).out_all_st;

fh = figure('color','w','OuterPosition',[100, 100, 1800, 800]);
n_col = 4; n_rows = 4; items = n_col * n_rows;
tl = tiledlayout(n_rows, n_col, 'TileSpacing', 'tight', 'Padding', 'loose');

diff = ll.tv - ll.ed;
[c, i_ll] = sort(diff);
% 1. Calculate the 98th percentile of the absolute differences
% This ensures 98% of your data falls within the color range
max_lim = prctile(abs(diff), 99.5);

% 2. Create the colormap
cmap_size = 256;
col_map = cbrewer2('RdBu', cmap_size);

% 3. Map each diff value to an index, CLAMPING values outside [-max_lim, max_lim]
% Using 'clamp' prevents indices < 1 or > 256 for extreme outliers
map_idx = round(interp1(linspace(-max_lim, max_lim, cmap_size), 1:cmap_size, diff, 'linear', 'extrap'));
map_idx = max(1, min(cmap_size, map_idx)); % Ensure bounds

% 4. Extract colors
bout_colors = col_map(map_idx, :);

selection = 'random';

if strcmp(selection, 'best')
    selected = i_ll(end - items + 2:end)';
elseif strcmp(selection, 'worst')
    selected = i_ll(1:items - 1)';
elseif strcmp(selection, 'random')
    % 1. Pick 16 random indices from the pool of available bout indices
    rand_indices = randperm(n_bouts, items - 1); 
    
    % 2. Extract the actual 'diff' values for these specific bouts
    rand_diffs = diff(rand_indices);
    
    % 3. Sort the random selection based on those 'diff' values
    [~, sort_idx] = sort(rand_diffs);
    selected = rand_indices(sort_idx); % Now they are ordered by diff

elseif strcmp(selection, 'selection')
    selected = [4 11];
    %selected = problematic_idx(end-7)';

elseif strcmp(selection, 'kl_based')
    [c_kl, i_kl] = sort(ll.js_div_partial);
    selected = i_kl(end - items + 2:end)';
    selected = i_kl(1:items - 1)';

end

for idx_selected_bout = selected %
 % Extract local row data
    freeze_row = freezes_ref(idx_selected_bout, :);
    freeze_ed_row = freezes_ref(idx_selected_bout, :);

    dur_s      = freeze_row.durations_s;
    is_censored = dur_s > censoring_val;
    dur_s(is_censored) = censoring_val + 1/60;

    % --- Model 1 (Time-Varying) ---
    cur_out_tv = out_all_tv(idx_selected_bout, :);
    [pdf_ddm_tv] = compute_pdf_tv_ddm(cur_out_tv, res_points);
    ll_tv(idx_selected_bout) = compute_likelihood(pdf_ddm_tv, is_censored, dur_s);

    % --- Model 1 (Stationary) ---
    cur_out_st = out_all_st(idx_selected_bout, :);
    [pdf_ddm_st] = compute_pdf_tv_ddm(cur_out_st, res_points);
    ll_st(idx_selected_bout) = compute_likelihood(pdf_ddm_st, is_censored, dur_s);

    % --- Model 1 (Theory) ---
    local_row = freeze_row;
    local_row.sm = local_row.sm_stat;
    nll_val_st = nll_fly_ddm_newer(data(1).est, local_row, res_points, ...
        strcat('model_', model_list{1}), 'iid', '', []);
    [~, f_ddm, fd_ddm] = nll_fly_ddm_newer(data(1).est, local_row, res_points, ...
        strcat('model_', model_list{1}), 'iid', 'p', []);
    ll_st_theory(idx_selected_bout) = -nll_val_st;

   % --- Model 2 (Extrema Detection) ---
    local_row_ed = freeze_ed_row;
    
    % FIX: Slice the data into a temporary variable first
    temp_soc_mot = data(2).extra_tv_n(idx_selected_bout, :); 
    
    % Define the struct entirely within the loop to ensure it is "temporary"
    signal_ed = struct(); 
    signal_ed.soc_mot_array = temp_soc_mot;

    [nll_val_ed] = nll_fly_ddm_newer(data(2).est, local_row_ed, res_points, ...
        strcat('model_', model_list{2}), 'iid', '', signal_ed);
    ll_ed(idx_selected_bout) = -nll_val_ed;

    % Ensure the second call also uses the sliced signal_ed
    [~, f_ed, fd_ed] = nll_fly_ddm_newer(data(2).est, local_row_ed, res_points, ...
        strcat('model_', model_list{2}), 'iid', 'p', signal_ed);

    nexttile
    hold on

    % 1. Create a mask for valid (post-truncation) data
    % Assumes pdf_ddm_tv.t is the time vector corresponding to the ddm vector
    valid_mask = pdf_ddm_tv.t >= res_points.truncation ;
    
    p_pdf = zeros(size(pdf_ddm_tv.ddm));
    q_pdf = zeros(size(pdf_ddm_tv.ddm));

    % 2. Extract only the valid PDF segments
    p_pdf(valid_mask) = pdf_ddm_tv.ddm(valid_mask);
    q_pdf(valid_mask) = pdf_ddm_st.ddm(valid_mask); % Extract local row data

    %plot(fd(1:end-1), f(1:end-1), 'LineWidth', 1.4, 'Color', col.stationary_sm)
    %stem(fd(end), f(end), 'Color', col.stationary_sm)
    plot(pdf_ddm_tv.t, p_pdf, 'LineWidth', 1.4, 'Color', col.timevarying_sm)
    stem(fd_ddm(end), pdf_ddm_tv.survival, 'Color', col.timevarying_sm)

%     plot(pdf_ddm_st.t, q_pdf, 'LineWidth', 1.4, 'Color', col.stationary_sm)
%     stem(fd_ddm(end), pdf_ddm_st.survival, 'Color', col.stationary_sm)

    plot(fd_ed(1:end-1), f_ed(1:end-1), 'LineWidth', 1.4, 'Color', col.extremadetection)
    stem(fd_ddm(end), f_ed(end), 'Color', col.extremadetection)

    ll_st_theory(idx_selected_bout) - ll_st(idx_selected_bout)
    plot(fd_ddm, f_ddm, 'k--');

    %trapz(pdf_ddm_tv.t(pdf_ddm_tv.t >= data(1).results.points.truncation), pdf_ddm_tv.ddm(pdf_ddm_tv.t >= data(1).results.points.truncation)) + pdf_ddm_tv.survival
    %trapz(pdf_ddm_st.t(pdf_ddm_st.t >= data(1).results.points.truncation), pdf_ddm_st.ddm(pdf_ddm_st.t >= data(1).results.points.truncation)) + pdf_ddm_st.survival

    ylabel('density')

    ax = gca;
    apply_generic(ax, 'xticks', [0 10]);

    if idx_selected_bout == selected(1)
        xlabel('Duration (s)')
        xticklabels([0 10])

    else
        xlabel('')
        xticklabels([])
    end

    % 1. Define your Y-limits first
    y_min = -1;
    y_max = 2.0;
    ylim(ax, [y_min, y_max]);
    y_range = y_max - y_min;

    % 2. Define Proportions (Adjust these percentages as needed)
    % Position imagesc at 10% from the bottom
    img_y_pos = y_min + (0.06 * y_range);
    % Position scatter at 5% below the bottom (outside axis)
    scatter_y_pos = y_min - (0.1 * y_range);
    % Position X-Label at 15% below the bottom
    xlabel_y_pos = y_min - (0.15 * y_range);

    % 3. Plotting
    create_vector_for_imagesc = nan(size(data(1).t));
    create_vector_for_imagesc(1:length(signal_ed.soc_mot_array)) = signal_ed.soc_mot_array';

    % Use img_y_pos for the Y-location of imagesc
    % Inside your for-loop...
    imagesc(data(1).t, img_y_pos, create_vector_for_imagesc, [0 1.5]);
    colormap(ax, col.vars.sm); % <--- SPECIFY 'ax' so it doesn't touch the global figure map

    colormap(col.vars.sm)
    xlim([0 11])

    ax.Clipping = 'on';
    set(ax, 'Layer', 'Top')
    ax.YAxis.Visible = 'off';

    % 4. Dynamic Positioning
    ax.XLabel.Position(2) = xlabel_y_pos;

    xline(dur_s, 'LineWidth', 1.2, 'LineStyle', '--', 'Color', 'k')

    % Use scatter_y_pos for the marker
    scatter(dur_s, scatter_y_pos, 125, 'filled', '^', ...
        'Clipping', 'off', 'MarkerFaceColor', bout_colors(idx_selected_bout, :), 'MarkerEdgeColor', 'k')

    
end


nexttile(items) 
hold on

% 1. Plot Histogram
[counts, centers] = histcounts(diff, -2:0.02:2);
max_c = max(counts);
start_pos = 0.2;
fill([centers(1), centers], [start_pos, (counts/max_c)*0.5 + start_pos, start_pos], ...
     'k', 'EdgeColor', 'none');

% --- MODIFIED SECTION START ---
% 2. Plot Colorbar Scale across the TOTAL range of the data
total_max = max(abs(diff));
x_range_total = linspace(-total_max, total_max, 512); 

% Plot the imagesc using the full range
imagesc(x_range_total, -0.4, x_range_total); 

% Pin the colormap limits to your 99.5 percentile
% This makes the colors "bleed" past the percentile limits
clim([-max_lim, max_lim]); 
colormap(gca, col_map); 

% 3. Add markers
for idx_selected_bout = selected
    val_diff = diff(idx_selected_bout); 
    scatter(val_diff, start_pos - 0.05, 40, 'filled', '^', ...
        'MarkerFaceColor', bout_colors(idx_selected_bout, :), ... 
        'MarkerEdgeColor', 'k', ...
        'Clipping', 'off');
end

% 4. Formatting
ax = gca;
% Set the visual limits slightly wider than the total data for breathing room
xlim([-total_max, total_max]); 
% --- MODIFIED SECTION END ---

ylim([-0.1, 1]);
xlabel('\Delta LL (tv - ed)', 'FontSize', 12);
set(gca, 'YTick', [], 'Box', 'off', 'Color', 'none');
apply_generic(ax, 'no_y', true, 'font_size', 20)


% 4. Save this separate legend
exporter(fh, paths_analysis, sprintf('examples_%s_%d.pdf', selection, selected(1)))
