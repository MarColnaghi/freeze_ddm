clearvars

% 3. Trial table (like y in mixture script)
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_tests/dddm', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
kde_estimates = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/bsl/kde_spontaneous', id_code, 'kde_estimates_bsl.mat'));
[~, idx] = unique(kde_estimates.Fkde, 'last');
extra.Fkde = kde_estimates.Fkde(idx); extra.xkde = kde_estimates.xkde(idx); extra.fkde = kde_estimates.fkde(idx); 

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20);
n_bouts = height(bouts_proc);

% Only save useful variables in the table
y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
y.smp = bouts_proc.avg_sm_freeze_norm;
predictors = y.Properties.VariableNames;

link_linear = @(x) x;     % log link for bound height
link_logistic = @(x) 1./(1 + exp(-x));     % log link for bound height

% 4. Model definition (ground truth via model struct)

model.mu = struct( ...
    'predictors', {{ struct('name', 'sm') }}, ...
    'ground_truth', 1.5, ...        % drift scale
    'link', link_linear );

model.lambda = struct( ...
    'predictors', {{ struct('name','intercept') }}, ...
    'ground_truth', 0.14, ...       % leak
    'link', link_linear );

model.theta = struct(...
    'predictors', {{ ...
    struct('name', 'fs'), ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', [0.24 0.8], ...
    'link', link_linear ...
    );

model.tndt = struct(...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', [0.09], ...
    'link', link_linear ...
    );

% 5. Evaluate model → trial-wise parameters
[gt, lbl] = get_ground_truth_vector(model);
gt_table  = array2table(gt, 'VariableNames', lbl);

ncomp_vars = evaluate_model(model, gt_table, y);

col = cmapper();

% Load the Motion Cache
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% 6. Fixed likelihood / PDE parameters (unchanged)
fixed.dt       = 1/60;
fixed.dx       = 0.01;
fixed.sigma_sq = 1.0;
fixed.x0       = 0.0;
fixed.T_trunc  = 0.3;
fixed.T_max    = 10.5;

% 7. Fine simulation parameters
sim.dt = 0.00005;

% 8. Time vectors
t_coarse = 0:fixed.dt:fixed.T_max;
t_fine   = 0:sim.dt:fixed.T_max;

% 9. Storage (matches second script style)
rt = table;
rt.tv = nan(n_bouts,1);

drifts_cell = cell(n_bouts,1);
trial_seed  = nan(n_bouts,1);

chunk_len = length(t_coarse);

% 10. Simulation loop (structurally identical to mixture script)

tic
for idx_bout = 1:n_bouts

    ons = bouts_proc.onsets(idx_bout);

    sum_motion = motion_cache(bouts_proc.fly(idx_bout));
    sm_raw{idx_bout} = sum_motion(ons:ons + chunk_len - 1) ./ 10;
    sm_chunk = sm_raw{idx_bout};

    % ---- Trial-wise parameters ----
    mu_st     = ncomp_vars.mu(idx_bout);
    lambda_i = ncomp_vars.lambda(idx_bout);
    theta_i  = ncomp_vars.theta(idx_bout);
    tndt_i  = ncomp_vars.tndt(idx_bout);

    % ---- Generate COARSE drift signal (model-visible) ----
    drifts_cell{idx_bout}  = sm_chunk;

    % ---- Apply drift scale ----
    sim_signal_coarse = sm_chunk * gt_table.mu_sm;

    % ---- Interpolate to FINE grid ----
    fine_signal = interp1(t_coarse(:), sim_signal_coarse(:), ...
                          t_fine(:), 'nearest');
    fine_signal(isnan(fine_signal)) = sim_signal_coarse(end);

    % ---- Simulation parameter vector (trial-wise) ----
    sim_params_vec = [ ...
        sim.dt, ...
        fixed.sigma_sq, ...
        fixed.x0, ...
        theta_i, ...
        lambda_i ];

    % ---- Unique RNG seed per trial ----
    unique_seed   = uint64(idx_bout * 1000 + randi([1 60]));
    trial_seed(idx_bout) = unique_seed;

    % ---- Run high-resolution DDM simulation ----
    res = sim_ddm_seeded(fine_signal, sim_params_vec, 1, unique_seed);

    % ---- Censoring ----
    if isnan(res)
        rt.tv(idx_bout) = fixed.T_max + fixed.dt + tndt_i;
    else
        rt.tv(idx_bout) = res + tndt_i;
    end

end
toc

% 11. Final output (same pattern as mixture script)
rt = [rt ncomp_vars];

% 12. Done
disp('Simulation complete.');

points.censoring = fixed.T_max;
points.truncation = fixed.T_trunc;

fh = figure('color', 'w', 'Position', [100, 100, 600, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
hold on
histogram(rt.tv, 1/120:1/20:points.censoring + fixed.dt * 3, 'EdgeColor', 'none', 'Normalization', 'pdf')
apply_generic(gca)
ylabel('Density'); xlabel('Duration (s)'); ax.YAxis.Visible = 'off';

% exporter(fh, paths, 'Durations.pdf')

extra.soc_mot_array = cell2mat(sm_raw)';

%  Now we added our vector column to the bouts table.
bouts_proc.durations_s = rt.tv;
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_1s_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.intercept = ones(height(y),1);

model_2_fit = 'sddmtv11';
model_results = run_fitting_newer(bouts_proc, points, model_2_fit, paths, 'export', true,  'ground_truth', gt_table, 'bads_display', true, 'pass_ndt', false, 'n_bads', 1, 'extra', extra);
plot_estimates('results', model_results, 'export', true, 'paths', paths)
