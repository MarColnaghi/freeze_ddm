clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% basis_projection — VISUAL AID
%
% Show how the raw person-period design matrix S (one column per lag τ, i.e.
% the sliding social-motion window built in design_matrix.m / the estimator)
% is PROJECTED onto a smooth raised-cosine basis B:
%
%       Xₖ  =  S · B
%     (r×K)   (r×nLag)(nLag×K)
%
% Left  : S  — the lagged history block, many highly-correlated lag columns.
% Below : B  — K raised-cosine bumps tiling the SAME lag axis (built exactly
%              as the estimator's Bsm = raised_cosine_basis(lag_fr, K, width)).
% Right : Xₖ — each row's history reduced to K smooth coefficients. This is the
%              design actually fed to the GLM: it collapses ~nLag collinear lags
%              onto a handful of well-separated, smooth predictors, so the
%              recovered kernel β(τ)=B·b is smooth and the fit is well-posed.
% ════════════════════════════════════════════════════════════════════════

% ── parameters (mirror design_matrix.m / hazard_kernel_sparse.m) ─────────────
fps             = 60;
kernel_past_s   = 4.0;    % seconds of causal history
kernel_future_s = 2;      % acausal control (shown here for context)
dt_frames       = 6;      % sample the hazard every 6 frames (0.1 s)
trunc_point     = 30;     % left truncation (0.5 s) = risk-set entry
entry_fr        = trunc_point;
contact_threshold = 70;   % px, collision censoring
grid_anchor     = 'offset';

nb_basis        = 12;     % # raised-cosine bumps (estimator uses nb_kernel+nb_acausal)
bump_width      = 1.25;   % bump half-width in units of centre spacing
mask_preonset   = true;  % true → zero out lags reaching before freeze onset (no pre-onset drive)
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

% ── smooth basis B and the projection Xk = S·B ───────────────────────────────
B  = raised_cosine_basis(lag_fr, nb_basis, bump_width);   % nLag × K, same call as the estimator
K  = size(B,2);
Xk = Sdisp * B;                                           % r × K  (NaN separator rows stay NaN)
fprintf('Projection: S is %d×%d, B is %d×%d, Xk is %d×%d  (%d lags -> %d smooth cols)\n', ...
    size(Sdisp,1), size(Sdisp,2), size(B,1), size(B,2), size(Xk,1), size(Xk,2), nLag, K);

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
set(axS,'XTick', -kernel_past_s:1:kernel_future_s, 'XTickLabel', []);
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
xline(0, 'k-', 'LineWidth', 2)

% ── B : the smooth basis, sharing the SAME lag axis, drawn under S ───────────
axB = nexttile(7, [1 2]); hold(axB,'on');
cmapB = cbrewer2('Spectral', K);           % kernel-basis colours (shared with the Xk x-ticks)
for j = 1:K
    plot(axB, t_ax, B(ord,j), '-', 'Color', cmapB(j,:), 'LineWidth', 2);
end
xlim(axB,[min(t_ax) max(t_ax)]); ylim(axB,[0 1.05]);
set(axB,'XTick', -kernel_past_s:1:kernel_future_s, 'YTick', [0 1]);
xlabel(axB,'lag relative to current frame (s)');
ylabel(axB,'weight');
apply_generic(axB, 'font_size', 24);
set(axB,'XTick', -kernel_past_s:1:kernel_future_s);
linkaxes([axS axB], 'x');
xline(0, 'k-', 'LineWidth', 2)

% ── Xk : the projected design (rows × K) ─────────────────────────────────────
% Order the columns by each bump's centre time so Xk reads past -> future,
% matching the S / B panels (raised_cosine_basis lays centres out by lag frame,
% i.e. future -> past, which is the reverse of the plotted t_rel axis).
c_lag  = linspace(min(lag_fr), max(lag_fr), K);   % bump centres in frames
c_trel = -c_lag / fps;                             % centre time rel. to un-freeze (s)
[~, bord] = sort(c_trel);                          % ascending t_rel: past (left) -> future (right)

axX = nexttile(3, [2 1]); hold(axX,'on');
imagesc(axX, 1:K, 1:nR, Xk(:,bord), 'AlphaData', ~isnan(Xk(:,bord)));
set(axX,'YDir','reverse','Color',[0.86 0.86 0.86]);
colormap(axX, redmap);
qx = max(1, prctile(abs(Xk(~isnan(Xk))), 99));
caxis(axX, [-qx qx]);
xlim(axX,[0.5 K+0.5]); ylim(axX,[0.5 nR+0.5]);
xlabel(axX,'basis (s)');
apply_generic(axX, 'font_size', 24);
xtickangle(45)
% one x-tick per column, coloured to match its basis bump (cmapB, in Xk column order)
lbl = cell(1,K);
for p = 1:K
    col = cmapB(bord(p),:);
    lbl{p} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%.1f', col(1), col(2), col(3), c_trel(bord(p)));
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
    fn = fullfile(outdir, sprintf('basis_projection_%s%s.png', grid_anchor, masktag));
    exportgraphics(fh, fn, 'Resolution', 200);
    [~, base] = fileparts(fn);
    exporter(fh, paths, [base '.pdf']);   % vector PDF -> paths.fig (repo convention)
    fprintf('\nSaved: %s\n       %s\n', fn, fullfile(paths.fig, [base '.pdf']));
end

% ════════════════════════════════════════════════════════════════════════
% local functions
% ════════════════════════════════════════════════════════════════════════
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
