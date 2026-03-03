
% Load the table with the computed log likelihoods
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};

paths_analysis = path_generator('folder', fullfile('momentary_integration','alternative_test', model_list{1}, run_list{1}));
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
tl = tiledlayout(n_rows, n_col, 'TileSpacing', 'tight', 'Padding', 'compact');

diff = ll.tv - ll.st;
[c, i] = sort(diff);
max_lim = max(abs(diff));
% 2. Create the colormap (e.g., 256 levels for smoothness)
cmap_size = 256;
col_map = cbrewer2('RdBu', cmap_size);

% 3. Map each diff value to an index in the colormap (1 to 256)
% This formula scales [-max_lim, max_lim] to [1, 256]
map_idx = round(interp1(linspace(-max_lim, max_lim, cmap_size), 1:cmap_size, diff));

% 4. Extract the specific colors for your bouts
bout_colors = col_map(map_idx, :);

selection = 'random';

if strcmp(selection, 'best')
    selected = i(end - items + 1:end)';
elseif strcmp(selection, 'worst')
    selected = i(1:items)';
elseif strcmp(selection, 'random')
    % 1. Pick 16 random indices from the pool of available bout indices
    rand_indices = randperm(n_bouts, items); 
    
    % 2. Extract the actual 'diff' values for these specific bouts
    rand_diffs = diff(rand_indices);
    
    % 3. Sort the random selection based on those 'diff' values
    [~, sort_idx] = sort(rand_diffs);
    selected = rand_indices(sort_idx); % Now they are ordered by diff
end

for idx_selected_bout = selected %

    % Current freeze row
    freeze_row = freezes_ref(idx_selected_bout, :);
    dur_s      = freeze_row.durations_s;
    is_censored = dur_s >= censoring_val;
    dur_s(is_censored) = censoring_val + 1/60;

    % Update external covariate for current bout
    ec.soc_mot_array = data.extra_tv(idx_selected_bout, :);

    % --- Model 1 (Integration) ---
    cur_out_tv = out_all_tv(idx_selected_bout, :);
    [pdf_ddm_tv] = compute_pdf_tv_ddm(cur_out_tv, data(1).results.points);
    ll_tv(idx_selected_bout) = compute_likelihood(pdf_ddm_tv, is_censored, dur_s);


    cur_out_st = out_all_st(idx_selected_bout, :);
    [pdf_ddm_st] = compute_pdf_tv_ddm(cur_out_st, data(1).results.points);
    ll_st(idx_selected_bout) = compute_likelihood(pdf_ddm_st, is_censored, dur_s);

    % --- Model 2 (Extrema Detection) ---
    %     %Using the prepared freeze_row and ec
    local_row = freeze_row;
    local_row.sm = local_row.sm_stat;

    nll_val_st = nll_fly_ddm_newer(data(1).est, local_row, res_points, ...
        strcat('model_', model_list{1}), 'iid', '', []);
    [n, f, fd] = nll_fly_ddm_newer(data(1).est, local_row, res_points, ...
        strcat('model_', model_list{1}), 'iid', 'p', []);
    ll_st_theory(idx_selected_bout) = -nll_val_st;

    nexttile
    hold on

    %plot(fd(1:end-1), f(1:end-1), 'LineWidth', 1.4, 'Color', col.stationary_sm)
    %stem(fd(end), f(end), 'Color', col.stationary_sm)
    plot(pdf_ddm_tv.t, pdf_ddm_tv.ddm, 'LineWidth', 1.4, 'Color', col.timevarying_sm)
    stem(fd(end), pdf_ddm_tv.survival, 'Color', col.timevarying_sm)

    plot(pdf_ddm_st.t, pdf_ddm_st.ddm, 'LineWidth', 1.4, 'Color', col.stationary_sm)
    stem(fd(end), pdf_ddm_st.survival, 'Color', col.stationary_sm)

    trapz(fd(1:end - 1), f(1:end - 1)) + f(end)
    trapz(pdf_ddm_tv.t(pdf_ddm_tv.t >= data(1).results.points.truncation), pdf_ddm_tv.ddm(pdf_ddm_tv.t >= data(1).results.points.truncation)) + pdf_ddm_tv.survival



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
    y_max = 5.0;
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
    create_vector_for_imagesc = nan(size(data.t));
    create_vector_for_imagesc(1:length(ec.soc_mot_array)) = ec.soc_mot_array';

    % Use img_y_pos for the Y-location of imagesc
    % Inside your for-loop...
    imagesc(data.t, img_y_pos, create_vector_for_imagesc, [0 1.5]);
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
        'Clipping', 'off', 'MarkerFaceColor', bout_colors(idx_selected_bout, :))

    
end


exporter(fh, paths_analysis, sprintf('examples_%s_%d.pdf', selection, selected(1)))

% --- Create a separate Legend/Scale Figure ---
f_leg = figure('Color', 'w', 'Position', [100, 100, 400, 200]);
tiledlayout(2, 1, 'TileSpacing', 'none')
nexttile(2)
ax_leg = gca; % Thin bar in the middle

% 1. Create a dummy imagesc to represent the RdBu scale
% We use a horizontal gradient from -max_lim to max_lim
x_range = linspace(-max_lim, max_lim, 256);
imagesc(x_range, [0 1], x_range); 
colormap(ax_leg, col_map); % Apply the RdBu map here
set(ax_leg, 'YTick', [], 'XAxisLocation', 'bottom');

% 2. Add the selected bouts as markers on this scale
hold on;
for k = 1:length(selected)
    idx_bout = selected(k);
    val_diff = diff(idx_bout);
    % Plot a marker at the exact difference value on the scale
    scatter(val_diff, -0.4, 65, 'filled', '^', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
end

% 3. Labels and formatting
xlabel('\Delta Log-Likelihood');
xlim([-max_lim, max_lim]);

apply_generic(ax_leg, 'no_y', true)

nexttile(1)
ax = gca;
histogram(c, -3:0.02:3, 'EdgeColor', 'none', 'FaceColor', 'k')
xlim([-max_lim, max_lim]);
apply_generic(ax, 'no_y', true, 'no_x', true)

% 4. Save this separate legend
exporter(f_leg, paths_analysis, sprintf('legend_%s_%d.pdf', selection, selected(1)));