% --- 1. Setup and Data Loading ---
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};

% Adjust these paths to match your local PhD directory structure if needed
paths_analysis = path_generator('folder', fullfile('momentary_integration','alternative_test', model_list{1}, run_list{1}));
paths_analysis.fig = fullfile(paths_analysis.fig, 'examples_motion');
col = cmapper('', 2);

% Load Likelihoods and Data
ll = struct2table(importdata(fullfile(paths_analysis.results, 'll_table.mat')));
data = alternative_test_collect_data(model_list, run_list);

% Load Motion Cache (from your analyse_sm logic)
% Assuming data(1).sim_params contains the path to motion_cache
motion_cache = importdata(fullfile(paths_analysis.cache_path, 'motion_cache.mat'));

% --- 2. Configuration ---
align_to = 'offset'; % 'onset' or 'offset'
trace_len = 600;    % Frames to show
selection_size = 500; % How many top/bottom bouts to include in the average
n_col = 4; n_rows = 4; items = n_col * n_rows;

% Remove the bouts that are censored?

freezes_ref = data(1).table;
n_bouts     = height(freezes_ref);
diff = ll.tv - ll.st;
[~, sort_idx] = sort(diff);

% --- 3. Figure Initialization ---
fh = figure('color','w','OuterPosition',[100, 100, 1800, 800]);
tl = tiledlayout(n_rows, n_col, 'TileSpacing', 'tight', 'Padding', 'compact');

% Define the two groups
bottom_idx = sort_idx(1:selection_size);                 % Best for Stationary (Red)
top_idx    = sort_idx(end-selection_size+1:end);         % Best for Time-Varying (Blue)

% --- 3. Process Traces (Matrix Allocation) ---
bottom_mat = nan(selection_size, trace_len);
top_mat    = nan(selection_size, trace_len);

if ~isfield(data(1).results.estimates_mean, 'tndt_intercept')
    data(1).results.estimates_mean.tndt_intercept = 0;
end

if strcmpi(align_to, 'onset'), x = 0:trace_len-1; else, x = -trace_len+1:0; end

% Fill Matrices
for i = 1:selection_size
    % Extract Bottom Group
    bottom_mat(i, :) = get_cropped_signal_local(motion_cache, data(1).table, bottom_idx(i), trace_len, data(1).results.estimates_mean, align_to);
    % Extract Top Group
    top_mat(i, :) = get_cropped_signal_local(motion_cache, data(1).table, top_idx(i), trace_len, data(1).results.estimates_mean, align_to);
end

% --- 4. Main Plot ---
fh = figure('color','w','Position',[100, 100, 600, 450]); 
hold on;

% Plot Individual Traces with high transparency (0.1)
plot(x, bottom_mat, 'Color', [hex2rgb(col.stationary_sm) 0.05], 'LineWidth', 0.5);
plot(x, top_mat, 'Color', [hex2rgb(col.timevarying_sm) 0.05], 'LineWidth', 0.5);

% Plot Grand Median Traces
plot(x, median(bottom_mat, 1, 'omitnan'), 'Color', col.stationary_sm, 'LineWidth', 2.5);
plot(x, median(top_mat, 1, 'omitnan'), 'Color', col.timevarying_sm, 'LineWidth', 2.5);

% Formatting
apply_generic(gca);
y_min = -3; y_max = 3;
ylim([y_min, y_max]);
if strcmpi(align_to, 'onset'), xlim([0 600]); else, xlim([-600 0]); end

% Labels
if strcmpi(align_to,'onset')
    xlabel({'Time (frames)', '[freeze start aligned]'});
    ylabel({'Soc. Motion Signal', '(start subtracted)'});
else
    xlabel({'Time (frames)', '[freeze break aligned]'});
    ylabel({'Soc. Motion Signal', '(endpoint subtracted)'});
end

% Add a simple legend for the grand averages
text(x(round(end*0.1)), 2.5, 'Best for Stationary', 'Color', col.stationary_sm, 'FontWeight', 'bold');
text(x(round(end*0.1)), 2.1, 'Best for Time-Varying', 'Color', col.timevarying_sm, 'FontWeight', 'bold');

% exporter(fh, paths_analysis, sprintf('grand_average_sm_%s.pdf', align_to));

% --- Helper for local execution ---
function cropped = get_cropped_signal_local(motion_cache, y, idx, trace_len, estimates_mean, align_to)
    signal = motion_cache(y.fly(idx)) ./ 10;
    cropped = nan(1, trace_len);
    onset = y.onsets(idx);
    dur_f = round(y.durations_s(idx) * 60 - estimates_mean.tndt_intercept * 60);
    if dur_f < 1; return; end
    
    end_idx = onset + dur_f - 1;
    if onset < 1 || end_idx > numel(signal); return; end
    
    seg = signal(onset:end_idx);
    n = min(numel(seg), trace_len);
    
    if strcmpi(align_to,'onset')
        cropped(1:n) = seg(1:n);
        if ~isnan(cropped(1)); cropped = cropped - cropped(1); end
    else
        cropped(end-n+1:end) = seg(end-n+1:end);
        if ~isnan(cropped(end)); cropped = cropped - cropped(end); end
    end
end