%% Script: Individual Freeze Likelihood Comparison
% Refactored to loop over bouts (freezes) first.

paths_exp = path_generator('folder', 'momentary_integration/alternative');

% --- 1. Configuration & Paths ---
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; 
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
col = cmapper('', 2);
ddm_params.kde_grid = 0:1/60:30;

% Model 1: Integration (DDM-based) | Model 2: Extrema Detection
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};
paths      = path_generator('folder', 'fitting_freezes/le');

% --- 2. Data Pre-loading & Preparation ---
data = struct();

for m = 1:length(model_list)
    res_folder = fullfile(paths.results, model_list{m}, run_list{m});
    
    extra   = importdata(fullfile(res_folder, 'extra.mat'));
    res     = importdata(fullfile(res_folder, sprintf('fit_results_%s.mat', model_list{m})));
    freezes_tv = importdata(fullfile(res_folder, 'surrogate.mat'));
    freezes_st = importdata(fullfile(res_folder, 'surrogate.mat'));

    % --- STEP A: Check equality BEFORE modification ---
    data(m).raw_table = freezes_tv; 
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
    freezes_tv.sm = extra.soc_mot_array;
    % Use the specific censoring point for this specific model run
    freezes_tv.durations(freezes_tv.durations_s >= res.points.censoring) = 631;
    
    data(m).extra_tv   = extra;
    data(m).results = res;
    data(m).table   = freezes_tv;
    data(m).func    = str2func(strcat('model_', model_list{m}));
    data(m).est     = table2array(res.estimates_mean(:, ~ismissing(res.estimates_mean)));
    data(m).out_all_tv = evaluate_model(data(m).func(), res.estimates_mean, freezes_tv);

    freezes_st.sm = freezes_st.sm * ones(1, 631);
    freezes_st.durations(freezes_st.durations_s >= res.points.censoring) = 631;

    data(m).extra_st   = freezes_st.sm;
    data(m).results = res;
    data(m).table   = freezes_st;
    data(m).func    = str2func(strcat('model_', model_list{m}));
    data(m).est     = table2array(res.estimates_mean(:, ~ismissing(res.estimates_mean)));
    data(m).out_all_st = evaluate_model(data(m).func(), res.estimates_mean, freezes_st);
end

% Consistency Check
% if ~isequaln(data(1).table, data(2).table)
%     error("The freeze tables from the two models do not match. Aborting.");
% else
%     fprintf('Freeze tables match. Starting bout-wise analysis...\n');
% end
%%

% --- 3. Main Loop: Iterate over Bouts ---
freezes_ref = data(1).table; % Reference table
n_bouts     = height(freezes_ref);
censoring_val = data(1).results.points.censoring;

ll_tv = nan(n_bouts, 1);
ll_st = nan(n_bouts, 1);
ll_st_2 = nan(n_bouts, 1);

for idx_bout = 1:n_bouts
    
    if mod(idx_bout, 50) == 0
        fprintf('Processing bout %d/%d\n', idx_bout, n_bouts);
        total_ll_tv = sum(ll_tv, 'omitnan');
        total_ll_st = sum(ll_st, 'omitnan');

        fprintf('Time-Varying Total LL: %.2f\n', total_ll_tv);
        fprintf('Stationary Total LL: %.2f\n', total_ll_st);

        fprintf('Percentage favouring: %.2f\n',  sum(ll_tv > ll_st) / idx_bout);
    end
    
    % Current freeze row
    freeze_row = freezes_ref(idx_bout, :);
    dur_s      = freeze_row.durations_s;
    is_censored = dur_s >= censoring_val;
    
    % Update external covariate for current bout
    ec.soc_mot_array = data(1).extra_tv.soc_mot_array(idx_bout, :);
   
    % --- Model 1 (Integration) ---
    out_tv = data(1).out_all_tv(idx_bout, :);
    [pdf_ddm_tv] = compute_pdf_tv_ddm(out_tv, data(1).results.points);
    ll_tv(idx_bout) = compute_likelihood(pdf_ddm_tv, is_censored, dur_s);
    
    out_st = data(1).out_all_st(idx_bout, :);
    [pdf_ddm_st] = compute_pdf_tv_ddm(out_st, data(1).results.points);
    ll_st(idx_bout) = compute_likelihood(pdf_ddm_st, is_censored, dur_s);

% 
%     rt_mex = nan(100000, 1);
%     for idx_sims = 1:100000
%         if rand > out_st.pmix
%             sim_params_vec = [1/60, 1, 0, out_tv.theta1, 0];
%             rt_mex(idx_sims) = sim_ddm_seeded(out_tv.mu1, sim_params_vec, 1, idx_sims);
%         else
%             sim_params_vec = [1/60, 1, 0, out_tv.theta2, 0];
%             rt_mex(idx_sims) = sim_ddm_seeded(out_tv.mu2, sim_params_vec, 1, idx_sims);
%         end
%     end
%     rt_mex(isnan(rt_mex)) = 11;
% 
%     % --- Model 2 (Extrema Detection) ---
%     % Using the prepared freeze_row and ec

    freeze_row.sm = freeze_row.sm(1); 
    [~, f, fd] = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
                                  strcat('model_', model_list{m}), 'iid', 'p', ec);
    nll_val_st = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
                                  strcat('model_', model_list{m}), 'iid', '', ec);
    ll_st_2(idx_bout) = -nll_val_st;
%

fh = figure('color','w','Position',[100, 100, 700, 300]);
plot([pdf_ddm_tv.t; 10.5 + 1/60], [pdf_ddm_tv.ddm; pdf_ddm_tv.survival], 'LineWidth', 2, 'DisplayName','PDE timevarying')
hold on
plot([pdf_ddm_tv.t; 10.5 + 1/60], [pdf_ddm_st.ddm; pdf_ddm_st.survival], 'LineWidth', 2, 'DisplayName','PDE stationary')
plot(fd, f, 'LineWidth', 2, 'DisplayName','theory stationary')
xline(dur_s)
apply_generic(gca)
legend('FontSize', 24,  'Box','off')
%plot(pdf_ddm_st.t, ec.soc_mot_array./10)

%     %histogram(rt_mex, [0:1/60:11], 'Normalization', 'pdf')
    drawnow
%     
%     drawnow
end

% --- 4. Final Output ---
total_ll_tv = sum(ll_tv, 'omitnan');
total_ll_st = sum(ll_st, 'omitnan');
total_ll_st_2 = sum(ll_st_2, 'omitnan');

fprintf('\nFinal Results:\n');
fprintf('Time-Varying Total LL: %.2f\n', total_ll_tv);
fprintf('Stationary Total LL: %.2f\n', total_ll_st);
fprintf('Stationary2 Total LL: %.2f\n', total_ll_st_2);

%%

fh = figure('color','w','Position',[100, 100, 1800, 800]);
n_col = 3; n_rows = 3; items = n_col * n_rows;
tl = tiledlayout(n_rows, n_col, 'TileSpacing', 'tight', 'Padding', 'tight');

diff = ll_tv - ll_st;
[c, i] = sort(diff);

selection = 'best';

if strcmp(selection, 'best')
    selected = i(end - items + 1:end);
elseif strcmp(selection, 'worst')
    selected = i(1:items);
end

for idx_selected_bout = [selected'] %[randi(n_bouts, items, 1)]' 

    % Current freeze row
    freeze_row = freezes_ref(idx_selected_bout, :);
    dur_s      = freeze_row.durations_s;
    is_censored = dur_s >= censoring_val;
    dur_s(is_censored) = censoring_val + 1/60;

    % Update external covariate for current bout
    ec.soc_mot_array = data(1).extra_tv.soc_mot_array(idx_selected_bout, :);

    % --- Model 1 (Integration) ---
    out_tv = data(1).out_all_tv(idx_selected_bout, :);
    [pdf_ddm_tv] = compute_pdf_tv_ddm(out_tv, data(1).results.points);
    ll_tv(idx_selected_bout) = compute_likelihood(pdf_ddm_tv, is_censored, dur_s);


    out_tv = data(1).out_all_st(idx_selected_bout, :);
    [pdf_ddm_st] = compute_pdf_tv_ddm(out_st, data(1).results.points);
    ll_st(idx_selected_bout) = compute_likelihood(pdf_ddm_st, is_censored, dur_s);

    % --- Model 2 (Extrema Detection) ---
%     %Using the prepared freeze_row and ec
%     [~, f, fd] = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
%                                   'model_dddm2', 'iid', 'p', ec);
%     nll_val_st = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
%                                   'model_dddm2', 'iid', '', ec);
%     ll_st(idx_bout) = -nll_val_st;

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
    
    xlabel('Duration (s)')
    ylabel('density')

    ax = gca;
    apply_generic(ax)

    create_vector_for_imagesc = nan(size(ddm_params.kde_grid));
    create_vector_for_imagesc(1:length(ec.soc_mot_array)) = ec.soc_mot_array;
    imagesc(ddm_params.kde_grid, -0.8, create_vector_for_imagesc, [0 1.2]);  % set XData = time_vector
    colormap(col.vars.sm)

    xlim([0 11])
    ylim([-1 6.0])
    xticks([0 10])

    ax.Clipping= 'on';
    set(ax ,'Layer', 'Top')
    ax.YAxis.Visible = 'off';
    ax.XLabel.Position(2) = -1.2;

    xline(dur_s, 'LineWidth', 1.2, 'LineStyle', '--', 'Color', 'k')
    scatter(dur_s, ax.YLim(1) - 0.15, 100, 'filled', '^', 'Clipping', 'off', 'MarkerFaceColor', 'k')

end

exporter(fh, paths_exp, sprintf('examples_%s_2.pdf', selection))

%%
fh = figure('color','w','Position',[100, 100, 1200, 300]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
ll = table();
ll.tv = ll_tv;
ll.st = ll_st_2;

total_ll_bar(ll, col)

nexttile(2, [1 2])
hold on
ylim([-3 3])
xlim([-50 height(ll) + 50])
ax = gca;

[sorted_deltall, idx_deltall_tv] = sort(ll.tv - ll.st);
idx_crossing = find(sorted_deltall > 0);

fill([ax.XLim(1) ax.XLim(1) idx_crossing(1) idx_crossing(1)], [ax.YLim(1) 0 0 ax.YLim(1)], hex2rgb(col.stationary_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)
fill([idx_crossing(1) idx_crossing(1) ax.XLim(2) ax.XLim(2)], [ax.YLim(2) 0 0 ax.YLim(2)], hex2rgb(col.timevarying_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)

xline(idx_crossing(1), 'Label', sprintf('%.2f%%', (sum(sorted_deltall > 0)./length(sorted_deltall)) * 100), 'FontSize', 20, 'LabelOrientation','horizontal')

bar(sorted_deltall, 'FaceColor', 'k', 'EdgeColor', 'none' )
xlabel('Sorted Freezes')
xticks([])
ylabel('$\log \!\big(\mathcal{L}_{\mathrm{tv}} / \mathcal{L}_{\mathrm{st}}\big)$', ...
    'Interpreter','latex')
apply_generic(ax)

exporter(fh, paths_exp, 'compare_total_ll.pdf')

function total_ll_bar(lls_output, col)

bh = bar([sum(lls_output.tv), sum(lls_output.st)], 'FaceColor', 'flat', 'EdgeColor', 'flat', 'LineWidth', 2);
bh.CData(1,:) = hex2rgb(col.timevarying_sm);
bh.CData(2,:) = hex2rgb(col.stationary_sm);

% if strcmp(gen_data, 'tv')
%         text(1, bh.YData(1), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)
% 
% elseif strcmp(gen_data, 'st')
%         text(2, bh.YData(2), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)
% 
% elseif strcmp(gen_data, 'ig')
%         text(3, bh.YData(3), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)
% 
% end

xticklabels({'acc', 'st'})
ylabel('Total log($\mathcal{L}$)', 'Interpreter', 'latex')
ax = gca;
apply_generic(ax)
ax.YAxis.Direction = 'reverse';
ylim([-15000 0]); %ax.YTickLabel(1) = {'worse'}; ax.YTickLabel(end) = {'better'};
text(2, ax.YLim(1) + 500, sprintf('$\\Delta log(\\mathcal{L}): %.2f$', sum(lls_output.tv) - sum(lls_output.st)), 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center')

end




%%%

% %% Script: Individual Freeze Likelihood Comparison
% % Refactored to loop over bouts (freezes) first.
% 
% paths_exp = path_generator('folder', 'momentary_integration/alternative');
% 
% % --- 1. Configuration & Paths ---
% threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; 
% id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
% col = cmapper('', 2);
% ddm_params.kde_grid = 0:1/60:30;
% 
% % Model 1: Integration (DDM-based) | Model 2: Extrema Detection
% model_list = {'dddim2'};
% run_list   = {'run01_280226_1053'};
% paths      = path_generator('folder', 'fitting_freezes/le');
% 
% % --- 2. Data Pre-loading & Preparation ---
% data = struct();
% 
% for m = 1:length(model_list)
%     res_folder = fullfile(paths.results, model_list{m}, run_list{m});
%     
%     extra   = importdata(fullfile(res_folder, 'extra.mat'));
%     res     = importdata(fullfile(res_folder, sprintf('fit_results_%s.mat', model_list{m})));
%     freezes = importdata(fullfile(res_folder, 'surrogate.mat'));
%     freezes_2 = importdata(fullfile(res_folder, 'surrogate.mat'));
% 
%     % --- STEP A: Check equality BEFORE modification ---
%     data(m).raw_table = freezes; 
% %     if m == 2
% %         if isequaln(data(1).raw_table, data(2).raw_table)
% %             fprintf('The two freeze tables match\n')
% %         else
% %             % If they don't match, you can use 'return' like your original script
% %             fprintf("The two freeze tables don't match \n")
% %             return 
% %         end
% %     end
% 
%     % --- STEP B: Now apply model-specific modifications ---
%     freezes.sm = extra.soc_mot_array;
%     % Use the specific censoring point for this specific model run
%     freezes.durations(freezes.durations_s >= res.points.censoring) = 631;
%     
%     data(m).extra   = extra;
%     data(m).results = res;
%     data(m).table   = freezes;
%     data(m).func    = str2func(strcat('model_', model_list{m}));
%     data(m).est     = table2array(res.estimates_mean(:, ~ismissing(res.estimates_mean)));
%     data(m).out_all = evaluate_model(data(m).func(), res.estimates_mean, freezes);
% 
%     freezes_2.sm = freezes_2.sm * ones(1, size(extra.soc_mot_array, 2));
%     freezes_2.durations(freezes_2.durations_s >= res.points.censoring) = 631;
% 
%     data(m).extra_2   = freezes_2.sm;
%     data(m).results = res;
%     data(m).table   = freezes;
%     data(m).func    = str2func(strcat('model_', model_list{m}));
%     data(m).est     = table2array(res.estimates_mean(:, ~ismissing(res.estimates_mean)));
%     data(m).out_all_2 = evaluate_model(data(m).func(), res.estimates_mean, freezes_2);
% end
% 
% % Consistency Check
% % if ~isequaln(data(1).table, data(2).table)
% %     error("The freeze tables from the two models do not match. Aborting.");
% % else
% %     fprintf('Freeze tables match. Starting bout-wise analysis...\n');
% % end
% %%
% % --- 3. Main Loop: Iterate over Bouts ---
% freezes_ref = data(1).table; % Reference table
% n_bouts     = height(freezes_ref);
% censoring_val = data(1).results.points.censoring;
% 
% ll_tv = nan(n_bouts, 1);
% ll_st = nan(n_bouts, 1);
% ll_st_2 = nan(n_bouts, 1);
% 
% for idx_bout = 1:n_bouts
%     
%     if mod(idx_bout, 100) == 0
%         fprintf('Processing bout %d/%d\n', idx_bout, n_bouts);
%         total_ll_tv = sum(ll_tv, 'omitnan');
%         total_ll_st = sum(ll_st, 'omitnan');
% 
%         fprintf('Time-Varying Total LL: %.2f\n', total_ll_tv);
%         fprintf('Stationary Total LL: %.2f\n', total_ll_st);
%     end
%     
%     % Current freeze row
%     freeze_row = freezes_ref(idx_bout, :);
%     dur_s      = freeze_row.durations_s;
%     is_censored = dur_s >= censoring_val;
%     
%     % Update external covariate for current bout
%     ec.soc_mot_array = data(1).extra.soc_mot_array(idx_bout, :);
%     
%     data(1).results.points.truncation = 0;
% 
%     % --- Model 1 (Integration) ---
%     out_f1 = data(1).out_all(idx_bout, :);
%     [pdf_ddm] = compute_pdf_tv_ddm(out_f1, data(1).results.points);
%     ll_tv(idx_bout) = compute_likelihood(pdf_ddm, is_censored, dur_s);
%     
%     out_f2 = data(1).out_all_2(idx_bout, :);
%     [pdf_ddm_2] = compute_pdf_tv_ddm(out_f2, data(1).results.points);
%     ll_st(idx_bout) = compute_likelihood(pdf_ddm_2, is_censored, dur_s);
% 
%     rt_mex = nan(100000, 1);
%     for idx_sims = 1:100000
%         if rand > out_f2.pmix
%             sim_params_vec = [1/60, 1, 0, out_f1.theta1, 0];
%             rt_mex(idx_sims) = sim_ddm_seeded(out_f1.mu1, sim_params_vec, 1, idx_sims);
%         else
%             sim_params_vec = [1/60, 1, 0, out_f1.theta1, 0];
%             rt_mex(idx_sims) = sim_ddm_seeded(out_f1.mu1, sim_params_vec, 1, idx_sims);
%         end
%     end
% 
%     rt_mex(isnan(rt_mex)) = 11;
% 
%     % --- Model 2 (Extrema Detection) ---
%     % Using the prepared freeze_row and ec
%     [~, f, fd] = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
%                                   strcat('model_', model_list{m}), 'iid', 'p', ec);
%     nll_val_st = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
%                                   strcat('model_', model_list{m}), 'iid', '', ec);
%     ll_st_2(idx_bout) = -nll_val_st;
% 
%     fh = figure('color','w','Position',[100, 100, 700, 300]);
%     plot(pdf_ddm.t, pdf_ddm.ddm)
%     hold on
%     plot(pdf_ddm.t, pdf_ddm_2.ddm)
%     plot(fd, f)
%     xline(dur_s)
%     histogram(rt_mex, [0:1/60:11], 'Normalization', 'pdf')
%     drawnow
%     
% %     drawnow
% end
% 
% % --- 4. Final Output ---
% total_ll_tv = sum(ll_tv, 'omitnan');
% total_ll_st = sum(ll_st, 'omitnan');
% 
% fprintf('\nFinal Results:\n');
% fprintf('Time-Varying Total LL: %.2f\n', total_ll_tv);
% fprintf('Stationary Total LL: %.2f\n', total_ll_st);
% 
% %%
% 
% fh = figure('color','w','Position',[100, 100, 1800, 800]);
% n_col = 3; n_rows = 3; items = n_col * n_rows;
% tl = tiledlayout(n_rows, n_col, 'TileSpacing', 'tight', 'Padding', 'tight');
% 
% diff = ll_tv - ll_st;
% [c, i] = sort(diff);
% 
% selection = 'worst';
% 
% if strcmp(selection, 'best')
%     selected = i(end - items + 1:end);
% elseif strcmp(selection, 'worst')
%     selected = i(1:items);
% end
% 
% for idx_selected_bout = [randi(n_bouts, items, 1)]' 
% 
%     % Current freeze row
%     freeze_row = freezes_ref(idx_selected_bout, :);
%     dur_s      = freeze_row.durations_s;
%     is_censored = dur_s >= censoring_val;
%     dur_s(is_censored) = censoring_val + 1/60;
% 
%     % Update external covariate for current bout
%     ec.soc_mot_array = data(1).extra.soc_mot_array(idx_selected_bout, :);
% 
%     % --- Model 1 (Integration) ---
%     out_f1 = data(1).out_all(idx_selected_bout, :);
%     [pdf_ddm] = compute_pdf_tv_ddm(out_f1, data(1).results.points);
%     ll_tv(idx_selected_bout) = compute_likelihood(pdf_ddm, is_censored, dur_s);
% 
%     % --- Model 2 (Extrema Detection) ---
%     % Using the prepared freeze_row and ec
%     [~, f, fd] = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
%                                   'model_dddm2', 'iid', 'p', ec);
%     nll_val_st = nll_fly_ddm_newer(data(1).est, freeze_row, data(1).results.points, ...
%                                   'model_dddm2', 'iid', '', ec);
%     ll_st(idx_bout) = -nll_val_st;
% 
%     nexttile
%     hold on
% 
%     plot(fd(1:end-1), f(1:end-1), 'LineWidth', 1.4, 'Color', col.stationary_sm)
%     stem(fd(end), f(end), 'Color', col.stationary_sm)
%     plot(pdf_ddm.t, pdf_ddm.ddm, 'LineWidth', 1.4, 'Color', col.timevarying_sm)
%     stem(fd(end), pdf_ddm.survival, 'Color', col.timevarying_sm)
% 
%     trapz(fd(1:end - 1), f(1:end - 1)) + f(end)
%     trapz(pdf_ddm.t(pdf_ddm.t >= data(1).results.points.truncation), pdf_ddm.ddm(pdf_ddm.t >= data(1).results.points.truncation)) + pdf_ddm.survival
%     
%     xlabel('Duration (s)')
%     ylabel('density')
% 
%     ax = gca;
%     apply_generic(ax)
% 
%     create_vector_for_imagesc = nan(size(ddm_params.kde_grid));
%     create_vector_for_imagesc(1:length(ec.soc_mot_array)) = ec.soc_mot_array;
%     imagesc(ddm_params.kde_grid, -0.501, create_vector_for_imagesc, [0 1.8]);  % set XData = time_vector
%     colormap(col.vars.sm)
% 
%     xlim([0 11])
%     ylim([-0.3 2.0])
%     xticks([0 10])
% 
%     ax.Clipping= 'on';
%     set(ax ,'Layer', 'Top')
%     ax.YAxis.Visible = 'off';
%     ax.XLabel.Position(2) = -0.45;
% 
%     xline(dur_s, 'LineWidth', 1.2, 'LineStyle', '--', 'Color', 'k')
%     scatter(dur_s, ax.YLim(1) - 0.15, 100, 'filled', '^', 'Clipping', 'off', 'MarkerFaceColor', 'k')
% 
% end
% 
% exporter(fh, paths_exp, sprintf('examples_%s.pdf', selection))
% 
% %%
% fh = figure('color','w','Position',[100, 100, 1200, 300]);
% tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact')
% nexttile
% ll = table();
% ll.tv = ll_tv;
% ll.st = ll_st;
% 
% total_ll_bar(ll, col)
% 
% nexttile(2, [1 2])
% hold on
% ylim([-3 3])
% xlim([-50 height(ll) + 50])
% ax = gca;
% 
% [sorted_deltall, idx_deltall_tv] = sort(ll.tv - ll.st);
% idx_crossing = find(sorted_deltall > 0);
% 
% fill([ax.XLim(1) ax.XLim(1) idx_crossing(1) idx_crossing(1)], [ax.YLim(1) 0 0 ax.YLim(1)], hex2rgb(col.stationary_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)
% fill([idx_crossing(1) idx_crossing(1) ax.XLim(2) ax.XLim(2)], [ax.YLim(2) 0 0 ax.YLim(2)], hex2rgb(col.timevarying_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)
% 
% xline(idx_crossing(1), 'Label', sprintf('%.2f%%', (sum(sorted_deltall > 0)./length(sorted_deltall)) * 100), 'FontSize', 20, 'LabelOrientation','horizontal')
% 
% bar(sorted_deltall, 'FaceColor', 'k', 'EdgeColor', 'none' )
% xlabel('Sorted Freezes')
% xticks([])
% ylabel('$\log \!\big(\mathcal{L}_{\mathrm{tv}} / \mathcal{L}_{\mathrm{st}}\big)$', ...
%     'Interpreter','latex')
% apply_generic(ax)
% 
% exporter(fh, paths_exp, 'compare_total_ll.pdf')
% 
% function total_ll_bar(lls_output, col)
% 
% bh = bar([sum(lls_output.tv), sum(lls_output.st)], 'FaceColor', 'flat', 'EdgeColor', 'flat', 'LineWidth', 2);
% bh.CData(1,:) = hex2rgb(col.timevarying_sm);
% bh.CData(2,:) = hex2rgb(col.stationary_sm);
% 
% % if strcmp(gen_data, 'tv')
% %         text(1, bh.YData(1), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)
% % 
% % elseif strcmp(gen_data, 'st')
% %         text(2, bh.YData(2), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)
% % 
% % elseif strcmp(gen_data, 'ig')
% %         text(3, bh.YData(3), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)
% % 
% % end
% 
% xticklabels({'acc', 'st'})
% ylabel('Total log($\mathcal{L}$)', 'Interpreter', 'latex')
% ax = gca;
% apply_generic(ax)
% ax.YAxis.Direction = 'reverse';
% ylim([-15000 0]); %ax.YTickLabel(1) = {'worse'}; ax.YTickLabel(end) = {'better'};
% text(2, ax.YLim(1) + 500, sprintf('$\\Delta log(\\mathcal{L}): %.2f$', sum(lls_output.tv) - sum(lls_output.st)), 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center')
% 
% end