%% Script: Individual Freeze Likelihood Comparison
paths_analysis = path_generator('folder', 'momentary_integration/alternative_test');

% --- 1. Configuration & Paths ---
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
col = cmapper('', 2);

% model_list = {'dddim2'};
% run_list   = {'run02_260325'};
model_list = {'dddm2'};
run_list   = {'run14_260325'};

paths_fit  = path_generator('folder', 'fitting_freezes/le_new');
data = alternative_test_collect_data(model_list, run_list);

% --- 2. Preparation ---
freezes_ref = data(1).table;
n_bouts     = height(freezes_ref);
censoring_val = data(1).results.points.censoring;
res_points = data(1).results.points; 

% Pre-allocate outputs
ll_tv = nan(n_bouts, 1);
ll_st = nan(n_bouts, 1);
ll_st_theory = nan(n_bouts, 1);
kl_tv_st = nan(n_bouts, 1);
kl_st_tv = nan(n_bouts, 1);
kl_tv_st_partial = nan(n_bouts, 1);
kl_st_tv_partial = nan(n_bouts, 1);
js_div = nan(n_bouts, 1);
js_div_partial = nan(n_bouts, 1);

out_all_tv = data(1).out_all_tv;
out_all_st = data(1).out_all_st;

% Ensure parallel pool is open
if isempty(gcp('nocreate')), parpool(6); end 

% --- 3. Parallel Checkpoint Setup ---
q = parallel.pool.DataQueue();
startTime = tic;

% Use an anonymous function to pass the arrays into the helper
afterEach(q, @(idx) helperUpdateProgress(ll_tv, ll_st, n_bouts, startTime));

fprintf('Starting parfor over %d bouts...\n', n_bouts);

% --- 4. Main Loop: parfor ---
for idx_bout = 1:n_bouts
    % Extract local row data
    freeze_row = freezes_ref(idx_bout, :);
    dur_s      = freeze_row.durations_s;
    is_censored = dur_s > censoring_val;

    % --- Model 1 (Time-Varying) ---
    cur_out_tv = out_all_tv(idx_bout, :);
    [pdf_ddm_tv] = compute_pdf_tv_ddm(cur_out_tv, res_points);
    ll_tv(idx_bout) = compute_likelihood(pdf_ddm_tv, is_censored, dur_s);

    % --- Model 1 (Stationary) ---
    cur_out_st = out_all_st(idx_bout, :);
    [pdf_ddm_st] = compute_pdf_tv_ddm(cur_out_st, res_points);
    ll_st(idx_bout) = compute_likelihood(pdf_ddm_st, is_censored, dur_s);

    % --- Model 2 (Extrema Detection) ---
    data(1).est(end)

    local_row = freeze_row;
    local_row.sm = local_row.sm_stat;
    nll_val_st = nll_fly_ddm_newer(data(1).est, local_row, res_points, ...
        strcat('model_', model_list{1}), 'iid', '', []);
    [~, f, fd] = nll_fly_ddm_newer(data(1).est, local_row, res_points, ...
        strcat('model_', model_list{1}), 'iid', 'p', []);
    ll_st_theory(idx_bout) = -nll_val_st;

    % 1. Create a mask for valid (post-truncation) data
    % Assumes pdf_ddm_tv.t is the time vector corresponding to the ddm vector
    valid_mask = pdf_ddm_tv.t >= res_points.truncation ;

    p_pdf = zeros(size(pdf_ddm_tv.ddm));
    q_pdf = zeros(size(pdf_ddm_tv.ddm));

    % 2. Extract only the valid PDF segments
    p_pdf(valid_mask) = pdf_ddm_tv.ddm(valid_mask);
    q_pdf(valid_mask) = pdf_ddm_st.ddm(valid_mask);

    % --- KL Divergence Calculation ---
    p_vec = [pdf_ddm_tv.ddm(:); pdf_ddm_tv.survival] + eps;
    q_vec = [pdf_ddm_st.ddm(:); pdf_ddm_st.survival] + eps;
    p_norm = sum(p_vec);
    q_norm = sum(q_vec);
    p_vec = p_vec / p_norm;
    q_vec = q_vec / q_norm;
    
    kl_tv_st(idx_bout) = sum(p_vec .* log2(p_vec ./ q_vec));
    kl_st_tv(idx_bout) = sum(q_vec .* log2(q_vec ./ p_vec));
    
    m_vec = 0.5 * (p_vec + q_vec);
    js_div(idx_bout) = 0.5 * sum(p_vec .* log2(p_vec ./ m_vec)) + ...
                       0.5 * sum(q_vec .* log2(q_vec ./ m_vec));
    
    % 2. Partial KL (To the break)
    time_v = [pdf_ddm_tv.t; 9000000];
    partial_mask = (time_v >= res_points.truncation) & (time_v <= dur_s);   
    p_part = p_vec(partial_mask);
    q_part = q_vec(partial_mask);
    
    % Normalize these segments to compare "shapes" during the freeze
    kl_tv_st_partial(idx_bout)  = sum(p_part .* log2(p_part ./ q_part));
    kl_st_tv_partial(idx_bout)  = sum(q_part .* log2(q_part ./ p_part));

    m_vec_partial = 0.5 * (p_part + q_part);

    js_div_partial(idx_bout) = 0.5 * sum(p_part .* log2(p_part ./ m_vec_partial)) + ...
        0.5 * sum(q_part .* log2(q_part ./ m_vec_partial));

    figure
    plot(fd,f )
    hold on
    plot(pdf_ddm_tv.t, p_pdf)
    plot(pdf_ddm_tv.t, q_pdf)
    plot(fd,f )
    plot(pdf_ddm_tv.t, cur_out_tv.mu1)
    % Update progress
    send(q, idx_bout);
end

% --- 5. Final Output & Saving ---
total_ll_tv = sum(ll_tv, 'omitnan');
total_ll_st = sum(ll_st, 'omitnan');
total_ll_st_theory = sum(ll_st_theory, 'omitnan');
fav_pct = (sum(ll_tv > ll_st, 'omitnan') / sum(~isnan(ll_tv))) * 100;

fprintf('\nFinal Results:\n');
fprintf('Time-Varying Total LL: %.2f\n', total_ll_tv);
fprintf('Stationary Total LL: %.2f\n', total_ll_st);
fprintf('Stationary2 Total LL: %.2f\n', total_ll_st_theory);
fprintf('Percentage favoring Time-Varying: %.2f%%\n', fav_pct);

% Folder versioning
results_base = fullfile(paths_analysis.results, model_list{1}, run_list{1});
version_idx = 1;
results_folder = results_base;
while exist(results_folder, 'dir')
    version_idx = version_idx + 1;
    results_folder = sprintf('%s_v%d', results_base, version_idx);
end
mkdir(results_folder);
cd(results_folder)

figure_base = fullfile(paths_analysis.fig, model_list{1}, run_list{1});
version_idx = 1;
figure_folder = figure_base;
while exist(results_folder, 'dir')
    version_idx = version_idx + 1;
    results_folder = sprintf('%s_v%d', figure_base, version_idx);
end

ll_table.st = ll_st;
ll_table.tv = ll_tv;
ll_table.st_theory = ll_st_theory;
ll_table.kl_tv_st = kl_tv_st;
ll_table.kl_st_tv = kl_st_tv;
ll_table.js_div = js_div;
ll_table.kl_tv_st_partial = kl_tv_st_partial;
ll_table.kl_st_tv_partial = kl_st_tv_partial;
ll_table.js_div_partial = js_div_partial;
save('ll_table.mat', 'll_table')

% --- Helper Function ---
function helperUpdateProgress(ll_tv, ll_st, n_bouts, startTime)
    persistent n_completed
    if isempty(n_completed)
        n_completed = 0;
    end
    n_completed = n_completed + 1;
    
    if mod(n_completed, 100) == 0
        elapsed = toc(startTime);
        % Note: Arrays are passed as they were at the time of the callback trigger
        current_fav_pct = (sum(ll_tv > ll_st, 'omitnan') / sum(~isnan(ll_tv))) * 100;
        total_diff = sum(ll_tv, 'omitnan') - sum(ll_st, 'omitnan');
        
        fprintf('\n--- Checkpoint: %d/%d Bouts --\n', n_completed, n_bouts);
        fprintf('Time: %.1f s | Favors TV: %.1f%% | LL Diff: %.2f\n', ...
            elapsed, current_fav_pct, total_diff);
    end
end