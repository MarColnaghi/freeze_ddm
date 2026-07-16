clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% basis_projection_exp — VISUAL AID  (companion to basis_projection.m)
%
% Same S · B = Xₖ projection, but B is the EXPONENTIAL BANK used by model (B)
% of hazard_kernel_sparse.m (exp_basis), not the raised-cosine basis:
%
%   • n_exp decaying exponentials exp(-τ/τ_k) on the CAUSAL side — all anchored
%     at "now" (t=0) and falling into the past at different time-constants τ_k
%     (log-spaced 50 ms → kernel_past_s). Each column is L2-normalised, so fast
%     decays are tall/narrow and slow decays are low/broad.
%   • nb_acausal raised-cosine bumps on the ACAUSAL (future) side — the leakage
%     control block (grey).
%
% Under the L1 (sparse) prior this bank turns the fit into TIMESCALE SELECTION:
% the retained τ_k are the leak estimate. The Xₖ columns are the atoms; each
% x-tick is coloured to its atom and labelled by τ (exp) or centre time (acausal).
% ════════════════════════════════════════════════════════════════════════

% ── parameters (mirror hazard_kernel_sparse.m) ───────────────────────────────
fps             = 60;
kernel_past_s   = 4.0;    % seconds of causal history
kernel_future_s = 2;      % acausal control (shown here for context)
dt_frames       = 6;      % sample the hazard every 6 frames (0.1 s)
trunc_point     = 30;     % left truncation (0.5 s) = risk-set entry
entry_fr        = trunc_point;
contact_threshold = 70;   % px, collision censoring
grid_anchor     = 'offset';

n_exp           = 8;      % # decaying-exponential atoms (bank; estimator uses ~10)
nb_acausal      = 4;      % # raised-cosine bumps on the acausal side (leakage control)
show_acausal    = false;  % false → leave the acausal block blank (axes unchanged)
bump_width      = 1.25;   % acausal bump half-width in units of centre spacing
mask_preonset   = true;   % true → zero out lags reaching before freeze onset (no pre-onset drive)
dur_lo = 123; dur_hi = 156;
save_out        = true;

% ── load + filter (identical to the real pipeline) ───────────────────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel/visual_aid','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);

bl = data_parser_new(bouts, 'type','immobility','period','loom','window','le','nloom',2:20,'min_dur',trunc_point);
bl = impose_contact_threshold(bl, 'threshold', contact_threshold, 'type','onlyfreeze');
bl.dur_frames = bl.ending_time;
bl = bl(bl.dur_frames > trunc_point, :);

% ── global z-scoring of sm ───────────────────────────────────────────────────
flies = unique(bl.fly); nn=0; s1=0; s2=0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn=nn+numel(v); s1=s1+sum(v); s2=s2+sum(v.^2);
end
sm_mu = s1/nn; sm_sd = sqrt(s2/nn - sm_mu^2);

% ── lag grid (same as the estimator) ─────────────────────────────────────────
back_fr = round(kernel_past_s*fps);  fwd_fr = round(kernel_future_s*fps);
lag_fr  = (-fwd_fr:back_fr)';  nLag = numel(lag_fr);  lag_s = lag_fr/fps;
t_rel   = -lag_s;                       % past < 0, now = 0, future > 0
caus_mask = lag_fr >= 0;                % causal (past) lags

% ── choose one uncensored + one censored bout (clean full window) ────────────
pick_unc = find_bout(bl, motion_cache, false, dur_lo, dur_hi, back_fr, fwd_fr, entry_fr);
pick_cen = find_bout(bl, motion_cache, true,  dur_lo, dur_hi, back_fr, fwd_fr, entry_fr);
if isempty(pick_unc) || isempty(pick_cen)
    error('Could not find both an uncensored and a censored bout in the duration range.');
end
sel = [pick_unc, pick_cen];

% ── build the at-risk rows S for the two bouts (same bookkeeping) ─────────────
Sdisp = []; boutRows = zeros(numel(sel),2); rlab = {};
onset_lag = []; offset_lag = [];   % per-row frames since onset / until offset (for the stairs)
for k = 1:numel(sel)
    b   = sel(k);
    sm  = motion_cache(bl.fly(b)); sm = sm(:); L = numel(sm);
    on  = bl.onsets(b); dur = bl.dur_frames(b);
    switch grid_anchor
        case 'entry',  base = entry_fr;
        case 'offset', base = dur;
        otherwise, error('grid_anchor must be ''entry'' or ''offset''.');
    end
    kk   = ceil((1-base)/dt_frames) : floor((dur-base)/dt_frames);
    allg = base + kk(:)*dt_frames;
    grid = allg(allg >= entry_fr);
    if strcmp(grid_anchor,'entry')
        if isempty(grid) || grid(end) < dur,     grid = [grid; dur];      end %#ok<AGROW>
    else
        if isempty(grid) || grid(1)  > entry_fr, grid = [entry_fr; grid]; end %#ok<AGROW>
    end
    Sk = []; tk = []; lok = []; ofk = [];
    for gi = 1:numel(grid)
        f = grid(gi); t_abs = on + f - 1;
        if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end
        s = (sm(t_abs - lag_fr) - sm_mu) / sm_sd;
        if any(isnan(s)), continue; end
        if mask_preonset, s(lag_fr > f-1) = 0; end   % drop pre-onset drive (0 = mean, neutral)
        Sk(end+1,:) = s'; tk(end+1,1) = f/fps; %#ok<AGROW>
        lok(end+1,1) = (f-1);    ofk(end+1,1) = (dur-f); %#ok<AGROW>  % frames since onset / until offset
    end
    r0 = size(Sdisp,1);
    Sdisp = [Sdisp; Sk]; %#ok<AGROW>
    onset_lag = [onset_lag; lok]; offset_lag = [offset_lag; ofk]; %#ok<AGROW>
    boutRows(k,:) = [r0+1, size(Sdisp,1)];
    for j = 1:numel(tk), rlab{end+1} = sprintf('%.2f', tk(j)); end %#ok<AGROW>
    if k < numel(sel)
        Sdisp = [Sdisp; nan(1,nLag)]; rlab{end+1} = ''; %#ok<AGROW>
        onset_lag = [onset_lag; NaN]; offset_lag = [offset_lag; NaN]; %#ok<AGROW>
    end
    fprintf('bout %d: fly=%d dur=%d (%.2fs) censored=%d -> %d rows\n', ...
        b, bl.fly(b), dur, dur/fps, bl.is_censored(b), size(Sk,1));
end
nR = size(Sdisp,1);

% ── exponential bank B and the projection Xk = S·B ───────────────────────────
tau_grid_s = logspace(log10(0.05), log10(kernel_past_s), n_exp);   % 50 ms → kernel_past_s
[B, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, nb_acausal, bump_width);
K  = size(B,2);
Xk = Sdisp * B;                                           % r × K  (NaN separator rows stay NaN)
% acausal bump centres (for the future-side tick labels)
ac_lag  = linspace(min(lag_fr(~caus_mask)), max(lag_fr(~caus_mask)), nb_acausal);
ac_trel = -ac_lag / fps;                                  % centre time (s), future > 0
fprintf('Projection: S is %d×%d, B is %d×%d (n_exp=%d + acausal=%d), Xk is %d×%d\n', ...
    size(Sdisp,1), size(Sdisp,2), size(B,1), size(B,2), n_exp, nb_acausal, size(Xk,1), size(Xk,2));

% ── colours: one per atom (exp = Spectral by τ; acausal = grey) ──────────────
cmapE   = cbrewer2('Spectral', n_exp);
colA    = [0.60 0.60 0.60];
atomcol = [cmapE; repmat(colA, nb_acausal, 1)];   % K×3

% ════════════════════════════════════════════════════════════════════════
% figure:   [ S ]  [ Xk ]
%           [ B ]
% ════════════════════════════════════════════════════════════════════════
fh = figure('Color','w','Position',[80 80 650 700]);
tl = tiledlayout(3, 3, 'TileSpacing','compact', 'Padding','compact');
redmap = cbrewer2('Reds', []);

% ── S : the raw design matrix (rows × lags) ──────────────────────────────────
axS = nexttile(1, [2 2]); hold(axS,'on');
[t_ax, ord] = sort(t_rel);                 % past on the LEFT, now/future on the RIGHT
imagesc(axS, t_ax, 1:nR, Sdisp(:,ord), 'AlphaData', ~isnan(Sdisp(:,ord)));
set(axS,'YDir','reverse','Color',[0.86 0.86 0.86]);
colormap(axS, redmap); caxis(axS,[-3 3]);
xlim(axS,[min(t_ax) max(t_ax)]); ylim(axS,[0.5 nR+0.5]);
ylabel(axS,'at-risk frame');
set(axS,'YTick',[]);
apply_generic(axS, 'font_size', 24);
set(axS,'XTick', -kernel_past_s:1:kernel_future_s, 'XTickLabel', []);

% freeze-onset (black) and freeze-offset (white) boundaries, as in design_matrix.m
by = []; bx = []; by2 = []; bx2 = [];
for r = 1:nR
    if isnan(onset_lag(r))
        by(end+1)=NaN; bx(end+1)=NaN; by2(end+1)=NaN; bx2(end+1)=NaN; continue; %#ok<AGROW>
    end
    by(end+1)=r;  bx(end+1)=-onset_lag(r)/fps;   %#ok<AGROW>  pre-onset is LEFT of this line
    by2(end+1)=r; bx2(end+1)=offset_lag(r)/fps;  %#ok<AGROW>  post-offset is RIGHT of this line
end
plot(axS, bx2, by2, 'w-.', 'LineWidth', 1);
plot(axS, bx,  by,  'k-.', 'LineWidth', 1);
xline(0, 'k-', 'LineWidth', 2);

% ── B : the exponential bank, sharing the SAME lag axis, drawn under S ────────
axB = nexttile(7, [1 2]); hold(axB,'on');
bmax = max(max(B(:, is_exp))) * 1.05;               % scale to the exponential atoms
for j = 1:K
    if ~show_acausal && ~is_exp(j), continue; end   % acausal left blank
    lw = 1.4 + 0.6*is_exp(j);
    plot(axB, t_ax, B(ord,j), '-', 'Color', atomcol(j,:), 'LineWidth', lw);
end
xlim(axB,[min(t_ax) max(t_ax)]); ylim(axB,[0 bmax]);
xlabel(axB,'lag relative to current frame (s)');
ylabel(axB,'weight');
apply_generic(axB, 'font_size', 24);
set(axB,'XTick', -kernel_past_s:1:kernel_future_s);
linkaxes([axS axB], 'x');
xline(0, 'k-', 'LineWidth', 2);

% ── Xk : the projected design (rows × K atoms; exp block | acausal block) ─────
% display order: reverse the exponential block so faster tau sits on the RIGHT;
% the acausal block stays at the far right (left blank when show_acausal is false)
perm    = [n_exp:-1:1, (n_exp+1):K];
Xk_disp = Xk(:, perm);
colp    = atomcol(perm, :);
isexp_p = is_exp(perm);
if ~show_acausal, Xk_disp(:, ~isexp_p) = NaN; end     % blank acausal columns (xlim unchanged)
axX = nexttile(3, [2 1]); hold(axX,'on');
imagesc(axX, 1:K, 1:nR, Xk_disp, 'AlphaData', ~isnan(Xk_disp));
set(axX,'YDir','reverse','Color',[0.86 0.86 0.86]);
colormap(axX, redmap);
qx = max(1, prctile(abs(Xk_disp(~isnan(Xk_disp))), 99));
caxis(axX, [-qx qx]);
xlim(axX,[0.5 K+0.5]); ylim(axX,[0.5 nR+0.5]);        % xlim spans all K atoms; acausal blank
xlabel(axX,'exponential timescale  \tau (s)');
apply_generic(axX, 'font_size', 24);
xtickangle(axX, 45);
if show_acausal
    plot(axX, (n_exp+0.5)*[1 1], [0.5 nR+0.5], 'k-', 'LineWidth', 0.8);   % exp | acausal divider
end
% one x-tick per atom, coloured to match its curve in B (acausal ticks left blank)
lbl = repmat({''}, 1, K);
for p = 1:K
    col = colp(p,:); orig = perm(p);
    if isexp_p(p)
        lbl{p} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%.2f', col(1),col(2),col(3), tau_grid(orig));
    elseif show_acausal
        a = orig - n_exp;
        lbl{p} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}+%.1f', col(1),col(2),col(3), ac_trel(a));
    end
end
set(axX,'XTick',1:K,'XTickLabel',lbl,'TickLabelInterpreter','tex','YTick',[]);
axX.XAxis.FontSize = 20;

cb = colorbar(axX); cb.Label.String = 'z-scored motion / projection';
cb.FontSize = 20;
cb.LineWidth = 2;

% ── save ─────────────────────────────────────────────────────────────────────
set([axS axB axX],'Layer','top');   % draw axis box / ticks on top of the image + overlays
axS.Toolbar.Visible='off'; axB.Toolbar.Visible='off'; axX.Toolbar.Visible='off';
if save_out
    outdir = fileparts(mfilename('fullpath'));
    masktag = ''; if mask_preonset, masktag = '_masked'; end
    fn = fullfile(outdir, sprintf('basis_projection_exp_%s%s.png', grid_anchor, masktag));
    exportgraphics(fh, fn, 'Resolution', 200);
    [~, base] = fileparts(fn);
    exporter(fh, paths, [base '.pdf']);   % vector PDF -> paths.fig (repo convention)
    fprintf('\nSaved: %s\n       %s\n', fn, fullfile(paths.fig, [base '.pdf']));
end

% ════════════════════════════════════════════════════════════════════════
% local functions
% ════════════════════════════════════════════════════════════════════════
function [B, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, nb_acausal, bump_width)
% Causal dictionary of decaying exponentials exp(-τ/τ_k) on τ≥0, plus a small
% raised-cosine block on the acausal (future) lags. Columns are L2-normalised so
% the L1 penalty is scale-fair across atoms (Mineault et al. 2009, Eq.11).
% Verbatim copy of exp_basis in hazard_kernel_sparse.m.
    lag_s = lag_fr(:) / fps;
    caus  = lag_s >= 0;
    ne    = numel(tau_grid_s);
    Be    = zeros(numel(lag_s), ne);
    for k = 1:ne
        col = zeros(numel(lag_s),1);
        col(caus) = exp(-lag_s(caus) / tau_grid_s(k));
        Be(:,k) = col / max(norm(col), eps);
    end
    Ba = zeros(numel(lag_s), nb_acausal);
    if any(~caus) && nb_acausal > 0
        Ba(~caus,:) = raised_cosine_basis(lag_fr(~caus), nb_acausal, bump_width);
        Ba = Ba ./ max(vecnorm(Ba), eps);
    end
    B        = [Be, Ba];
    tau_grid = tau_grid_s;
    is_exp   = [true(1,ne), false(1,nb_acausal)];
end

function idx = find_bout(bl, motion_cache, want_cens, dlo, dhi, back_fr, fwd_fr, entry_fr)
    idx = [];
    for b = 1:height(bl)
        if logical(bl.is_censored(b)) ~= want_cens, continue; end
        dur = bl.dur_frames(b); if dur < dlo || dur > dhi, continue; end
        on = bl.onsets(b); sm = motion_cache(bl.fly(b)); sm = sm(:); L = numel(sm);
        t0 = on + entry_fr - 1;
        if (t0 - back_fr) < 1 || (on + dur - 1 + fwd_fr) > L, continue; end
        w = sm((t0-back_fr):(on+dur-1+fwd_fr));
        if any(isnan(w)), continue; end
        idx = b; return
    end
end
