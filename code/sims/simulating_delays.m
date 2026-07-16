%% ════════════════════════════════════════════════════════════════════════
%  SIMULATING DELAYS — what each delay parameter actually does
%  ------------------------------------------------------------------------
%  The model has three delays that are easy to confuse because they all push
%  reaction times later. They act through completely different mechanisms:
%
%    sensory_delay  -> shifts the INPUT SIGNAL (evidence arrives late; the
%                      pre-onset interval is zero-padded). Changes the drive.
%    delayed_start  -> leaves the signal untouched but delays WHEN THE
%                      ACCUMULATOR STARTS integrating. Changes the path.
%    ndt            -> touches neither signal nor path; a pure POST-HOC shift
%                      of the reported RT. Changes only the readout.
%
%  The figure separates them: one column per delay, and three rows showing
%  where each one becomes visible.
%    Row 1 (input drift)   -> only sensory_delay moves it.
%    Row 2 (trajectories)  -> only sensory_delay and delayed_start move it.
%    Row 3 (RT distribution) -> all three move it.
%  Trajectories use the SAME SEEDS across every condition, so any difference
%  is the manipulation, not the noise.
%
%  Driven by a real freeze's social-motion trace.
% ════════════════════════════════════════════════════════════════════════
clearvars; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'simulators'));
addpath(fullfile(this_dir, 'simulators', 'cpp_mex_code'));

% ── Model regime ─────────────────────────────────────────────────────────
dt     = 1/60;          % s / frame (native social-motion rate)
sigma  = 1;           % diffusion SD
theta  = 2.7;           % bound -> median RT ~1.1 s on this trial (delays visible)
beta   = 4.2;           % drift gain (rate per unit social motion)
lambda = 0;             % leak (0 = perfect integrator)

idx_trial   = 5422;412;6000;%5350;%1223;     % a real freeze with a healthy motion trace
N           = 2000000;   % trials for the RT distributions
seeds_traj  = [74 211 300];   % example trajectories (SAME seeds in every condition)

delay_vals = [0 0.4 0.8];                              % s
delay_types = {'sensory_delay', 'delayed_start', 'ndt'};
delay_names = {'sensory delay', 'delayed start of accumulation', 'non-decision time'};
cols = [0.20 0.20 0.20; 0.90 0.50 0.13; 0.20 0.45 0.75];

% ── Real social-motion drive (one freeze) ────────────────────────────────
paths      = path_generator('folder', 'sims/sanity_checks');
bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type','immobility', 'period','loom', ...
                             'window','le', 'nloom', 2:20);
bouts_proc.smp       = ones(height(bouts_proc), 1);
bouts_proc.intercept = ones(height(bouts_proc), 1);
sm_signal  = extract_sm_from_bouts(bouts_proc, 'type','onsets', ...
                             'output_type','mat', 'window', [0 630], 'norm_factor', 10);

sm = sm_signal(idx_trial, :).';
sm = fillmissing(sm, 'previous');      % guard against motion gaps
sm(isnan(sm)) = 0;
n_fr  = numel(sm);
t_sig = (0:n_fr-1)' * dt;              % freeze-relative time, t=0 at freeze onset

base = struct('dt', dt, 'theta', theta, 'sigma', sigma, 'leak', lambda, ...
              'gain', beta, 'bias', 0, 'backend', 'mex', 'seed', 42);

xmax = 10.5;   % s, plotting window

fprintf('Trial %d | theta=%.2f beta=%.2f | %d frames (%.1f s)\n', ...
        idx_trial, theta, beta, n_fr, n_fr*dt);

%% ── Figure: 3 delays x 3 views ──────────────────────────────────────────
figure('Color','w','Position',[40 40 1500 950]);
tiledlayout(3, 3, 'TileSpacing','compact', 'Padding','compact');

RT = cell(numel(delay_types), numel(delay_vals));

sig_ylim  = [0 12];
traj_ylim = [-2 theta+0.5];
pdf_ylim  = [0 2];
% subtitle() is placed by the layout itself, so the note can never drift onto a
% neighbouring panel the way a manually positioned text() can.
note = @(ax, str) subtitle(ax, str, 'FontSize', 12, 'FontAngle','italic', ...
                           'Color', [0.35 0.35 0.35]);

% ---------- Row 1: the input drift the accumulator actually sees ----------
for c = 1:numel(delay_types)
    ax = nexttile(c); hold(ax, 'on')
    for v = 1:numel(delay_vals)
        p = base; p.(delay_types{c}) = delay_vals(v);
        drift = signal_to_drift(sm, p);
        plot(ax, t_sig, drift, 'Color', [cols(v,:) 0.9], 'LineWidth', 1.6, ...
             'DisplayName', sprintf('%.2f s', delay_vals(v)));
        % mark where accumulation begins (only delayed_start moves it)
        if strcmp(delay_types{c}, 'delayed_start')
            xline(ax, delay_vals(v), '--', 'Color', cols(v,:), 'LineWidth', 1.4, ...
                  'HandleVisibility','off');
        end
    end
    xlabel(ax, 'time from freeze onset (s)');
    if c == 1, ylabel(ax, 'drift  \beta\cdot sm(t)'); end
    title(ax, delay_names{c}, 'FontWeight','normal');
    switch c
        case 1
            legend(ax, 'Location','northeast','Box','off','FontSize',12);
            note(ax, 'signal shifts right + zero-pads');
        case 2, note(ax, 'signal identical; dashed = accumulation start');
        case 3, note(ax, 'signal identical (curves overlap)');
    end
    apply_generic(ax, 'xlim', [0 xmax], 'ylim', sig_ylim, 'font_size', 15, 'line_width', 1.5);
end

% ---------- Row 2: example trajectories (same seeds everywhere) ----------
for c = 1:numel(delay_types)
    ax = nexttile(3 + c); hold(ax, 'on')
    for v = 1:numel(delay_vals)
        p = base; p.(delay_types{c}) = delay_vals(v);
        drift  = signal_to_drift(sm, p);
        dstart = 0; if strcmp(delay_types{c}, 'delayed_start'), dstart = delay_vals(v); end
        ndt_c  = 0; if strcmp(delay_types{c}, 'ndt'),           ndt_c  = delay_vals(v); end
        n_skip = round(dstart / dt);
        drift_sim = drift(n_skip+1 : end);

        for s = 1:numel(seeds_traj)
            [k, traj] = sim_leaky_accumulator(drift_sim, theta, sigma, lambda, ...
                                              dt, seeds_traj(s), 0);
            t_abs = dstart + (1:numel(traj))' * dt;     % back to freeze-relative time
            plot(ax, t_abs, traj, '-', 'Color', [cols(v,:) 0.75], 'LineWidth', 1.4);
            if ~isnan(k)
                t_cross = dstart + k*dt;
                plot(ax, t_cross, theta, 'o', 'MarkerSize', 6, 'MarkerFaceColor', cols(v,:), ...
                     'MarkerEdgeColor','none');
                if ndt_c > 0   % ndt: path is identical, only the readout moves
                    plot(ax, [t_cross t_cross+ndt_c], [theta theta], ':', ...
                         'Color', cols(v,:), 'LineWidth', 1.6);
                    plot(ax, t_cross+ndt_c, theta, 'v', 'MarkerSize', 7, ...
                         'MarkerFaceColor', cols(v,:), 'MarkerEdgeColor','none');
                end
            end
        end
    end
    yline(ax, theta, 'k:', '\theta', 'LineWidth', 1.4, 'FontSize', 12);
    xlabel(ax, 'time from freeze onset (s)');
    if c == 1, ylabel(ax, 'DV'); end
    switch c
        case 1, note(ax, 'same noise, later evidence');
        case 2, note(ax, 'path starts later');
        case 3, note(ax, 'paths identical; \nabla = reported RT');
    end
    apply_generic(ax, 'xlim', [0 xmax], 'ylim', traj_ylim, 'font_size', 15, 'line_width', 1.5);
end

% ---------- Row 3: full RT distributions ----------
% bins aligned to the dt grid (RTs are exact multiples of dt) to avoid aliasing
edges = -dt/2 : 3*dt : xmax;
for c = 1:numel(delay_types)
    ax = nexttile(6 + c); hold(ax, 'on')
    for v = 1:numel(delay_vals)
        p = base; p.(delay_types{c}) = delay_vals(v);
        rt = simulate_freeze_ddm(sm, p, N);
        RT{c, v} = rt;
        r = rt(~isnan(rt));
        histogram(ax, r, edges, 'Normalization','pdf', 'EdgeColor','none', ...
                  'FaceColor', cols(v,:), 'FaceAlpha', 0.35, ...
                  'DisplayName', sprintf('%.2f s  (med %.2f)', delay_vals(v), median(r)));
        xline(ax, median(r), '--', 'Color', cols(v,:), 'LineWidth', 1.6, ...
              'HandleVisibility','off');
    end
    xlabel(ax, 'reaction time (s)');
    if c == 1, ylabel(ax, 'pdf'); end
    legend(ax, 'Location','northeast','Box','off','FontSize',11);
    note(ax, 'all three shift the RT');
    apply_generic(ax, 'xlim', [0 xmax], 'ylim', pdf_ylim, 'font_size', 15, 'line_width', 1.5);
end

%% ── Console summary: median shift per delay ─────────────────────────────
fprintf('\nMedian RT (s) and shift vs the 0 s baseline:\n');
fprintf('%-32s', 'delay');
fprintf('%10.2f s', delay_vals); fprintf('\n');
for c = 1:numel(delay_types)
    meds = cellfun(@(r) median(r(~isnan(r))), RT(c,:));
    fprintf('%-32s', delay_names{c});
    fprintf('%12.3f', meds); fprintf('\n');
    fprintf('%-32s', '   shift vs baseline');
    fprintf('%12.3f', meds - meds(1)); fprintf('\n');
end
fprintf(['\nAll three shift the RT by ~the delay, but only sensory_delay changes the\n' ...
         'drive (row 1) and only sensory_delay/delayed_start change the path (row 2).\n']);
