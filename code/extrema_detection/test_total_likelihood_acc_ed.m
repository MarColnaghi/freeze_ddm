%% Script: Individual Freeze Likelihood Comparison
% Refactored to loop over bouts (freezes) first.

% --- 1. Configuration & Paths ---
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; 
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
col = cmapper('', 2);

% Model 1: Integration (DDM-based) | Model 2: Extrema Detection
model_list = {'dddm2', 'ded2'};
run_list   = {'run21', 'run15'};
paths      = path_generator('folder', 'fitting_freezes/le');

% --- 2. Data Pre-loading & Preparation ---
data = struct();

for m = 1:length(model_list)
    res_folder = fullfile(paths.results, model_list{m}, run_list{m});
    
    extra   = importdata(fullfile(res_folder, 'extra.mat'));
    res     = importdata(fullfile(res_folder, sprintf('fit_results_%s.mat', model_list{m})));
    freezes = importdata(fullfile(res_folder, 'surrogate.mat'));
    
    % --- STEP A: Check equality BEFORE modification ---
    data(m).raw_table = freezes; 
%     if m == 2
%         if isequaln(data(1).raw_table, data(2).raw_table)
%             fprintf('The two freeze tables match\n')
%         else
%             % If they don't match, you can use 'return' like your original script
%             fprintf("The two freeze tables don't match \n")
%             return 
%         end
%     end

    % --- STEP B: Now apply model-specific modifications ---
    freezes.sm = extra.soc_mot_array;
    % Use the specific censoring point for this specific model run
    freezes.durations(freezes.durations_s >= res.points.censoring) = 631;
    
    data(m).extra   = extra;
    data(m).results = res;
    data(m).table   = freezes;
    data(m).func    = str2func(strcat('model_', model_list{m}));
    data(m).est     = table2array(res.estimates_mean(:, ~ismissing(res.estimates_mean)));
    data(m).out_all = evaluate_model(data(m).func(), res.estimates_mean, freezes);
end

% Consistency Check
% if ~isequaln(data(1).table, data(2).table)
%     error("The freeze tables from the two models do not match. Aborting.");
% else
%     fprintf('Freeze tables match. Starting bout-wise analysis...\n');
% end

% --- 3. Main Loop: Iterate over Bouts ---
freezes_ref = data(1).table; % Reference table
n_bouts     = height(freezes_ref);
censoring_val = data(1).results.points.censoring;

ll_integration = nan(n_bouts, 1);
ll_extrema     = nan(n_bouts, 1);

for idx_bout = 1:n_bouts
    if mod(idx_bout, 100) == 0, fprintf('Processing bout %d/%d\n', idx_bout, n_bouts); end
    
    % Current freeze row
    freeze_row = freezes_ref(idx_bout, :);
    dur_s      = freeze_row.durations_s;
    is_censored = dur_s >= censoring_val;
    
    % Update external covariate for current bout
    ec.soc_mot_array = data(1).extra.soc_mot_array(idx_bout, :);
    
    % --- Model 1 (Integration) ---
    out_f1 = data(1).out_all(idx_bout, :);
    [pdf_ddm] = compute_pdf_tv_ddm(out_f1, data(1).results.points);
    ll_integration(idx_bout) = compute_likelihood(pdf_ddm, is_censored, dur_s);
    
    % --- Model 2 (Extrema Detection) ---
    % Using the prepared freeze_row and ec
    [~, f, fd] = nll_fly_ddm_newer(data(2).est, freeze_row, data(2).results.points, ...
                                  'model_ded2', 'iid', 'p', ec);
    nll_val_ed = nll_fly_ddm_newer(data(2).est, freeze_row, data(2).results.points, ...
                                  'model_ded2', 'iid', '', ec);
    ll_extrema(idx_bout) = -nll_val_ed;

    fh = figure('color','w','Position',[100, 100, 750, 700]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
    nexttile

    hold on
    plot(pdf_ddm.t, pdf_ddm.ddm, 'Color', col.timevarying_sm, 'LineWidth', 2)
    stem(points.censoring, pdf_ddm.survival, 'Color', col.timevarying_sm, 'LineWidth', 2)
    plot(fd, f, 'Color', col.extremadetection, 'LineWidth', 2)
    apply_generic(gca, 'ylim', [0 0.025])
    xline(points.truncation, 'LineWidth', 3)
    xline(dur_s)

    text(5, 0.02, sprintf('difference: %.2f', ll_integration(idx_bout) - ll_extrema(idx_bout)))
    nexttile
    plot(pdf_ddm.t, ec.soc_mot_array,'Color', col.vars.sm(2, :) )
    xline(dur_s)
    apply_generic(gca)
    xline(points.truncation, 'LineWidth', 3)

    sum(pdf_ddm.ddm(pdf_ddm.t >= points.truncation)) + pdf_ddm.survival
    sum(f)
    hold off

end

% --- 4. Final Output ---
total_ll_int = sum(ll_integration, 'omitnan');
total_ll_ext = sum(ll_extrema, 'omitnan');

fprintf('\nFinal Results:\n');
fprintf('Integration Total LL: %.2f\n', total_ll_int);
fprintf('Extrema Det. Total LL: %.2f\n', total_ll_ext);