clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% design_matrix — VISUAL AID
%
% Show the person-period hazard-GLM design matrix for a FEW real freeze bouts
% (one that un-freezes, one censored by collision), built with the EXACT same
% pipeline / parameters as the real analysis (hazard_kernel_sparse.m):
%   • same loading + filtering (imm2_mob2_pc4, contact threshold, truncation)
%   • same lag grid (kernel_past_s of causal history, sampled at fps)
%   • same person-period sampling (entry_fr left truncation, every dt_frames)
%   • same global z-scoring of sm
% The main panel is the raw lagged-history block S (one row per at-risk frame,
% one column per lag τ) — i.e. the sliding social-motion window. A companion
% strip shows the outcome y (un-freeze event vs. censored). The dashed line marks
% the FREEZE-ONSET boundary: cells to its right are pre-onset (before the freeze
% began). Flip mask_preonset to zero those out (the within-freeze-only variant).
% ════════════════════════════════════════════════════════════════════════

% ── parameters (mirror hazard_kernel_sparse.m) ───────────────────────────────
fps             = 60;
kernel_past_s   = 4.0;    % seconds of causal history  (240 lags)
kernel_future_s = 2;      % acausal control (0 in the sparse driver)
dt_frames       = 6;      % sample the hazard every 6 frames (0.1 s)
trunc_point     = 30;     % left truncation (0.5 s); also the risk-set entry time
entry_fr        = trunc_point;
contact_threshold = 70;   % px, collision censoring

mask_preonset   = true;  % false = EXACTLY the real analysis; true = zero pre-onset lags
grid_anchor     = 'offset';  % 'entry'  = bins from risk-set entry (short bin lands on the EVENT)
                             % 'offset' = bins from the un-freeze/censor (regular event bin; short bin at entry)
dur_lo = 123; dur_hi = 156; % pick bouts in this duration range (keeps the image small)
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
fprintf('Bouts after filtering: %d (%.1f%% censored)\n', height(bl), 100*mean(bl.is_censored));

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

% ── choose one UNCENSORED and one CENSORED bout (clean full window) ───────────
pick_unc = find_bout(bl, motion_cache, false, dur_lo, dur_hi, back_fr, fwd_fr, entry_fr);
pick_cen = find_bout(bl, motion_cache, true,  dur_lo, dur_hi, back_fr, fwd_fr, entry_fr);
if isempty(pick_unc) || isempty(pick_cen)
    error('Could not find both an uncensored and a censored bout in the duration range.');
end
sel = [pick_unc, pick_cen];

% ── build the person-period rows for the two bouts (exact real bookkeeping) ───
rows_S = {}; rows_y = {}; rows_tinf = {}; rows_lagonset = {}; rows_offlag = {}; boutmeta = {};
rows_Str = {}; rows_ttr = {}; rows_lotr = {}; rows_oftr = {};   % pre-entry (truncated) rows
for k = 1:numel(sel)
    b   = sel(k);
    sm  = motion_cache(bl.fly(b)); sm = sm(:); L = numel(sm);
    on  = bl.onsets(b); dur = bl.dur_frames(b);
    % sampling grid over [entry_fr, dur], phase-locked to the chosen anchor
    switch grid_anchor
        case 'entry',  base = entry_fr;   % bins aligned to risk-set ENTRY -> short bin lands on the event
        case 'offset', base = dur;        % bins aligned to the OFFSET/event -> short bin lands at entry
        otherwise, error('grid_anchor must be ''entry'' or ''offset''.');
    end
    kk   = ceil((1-base)/dt_frames) : floor((dur-base)/dt_frames);
    allg = base + kk(:)*dt_frames;                       % phase-locked frames within [1, dur]
    grid = allg(allg >= entry_fr);
    if strcmp(grid_anchor,'entry')
        if isempty(grid) || grid(end) < dur,     grid = [grid; dur];      end   % exact OFFSET frame (event bin)
    else
        if isempty(grid) || grid(1)  > entry_fr, grid = [entry_fr; grid]; end   % exact ENTRY frame (entry bin)
    end
    Sk = []; yk = []; tk = []; lok = []; ofk = [];
    for gi = 1:numel(grid)
        f = grid(gi); t_abs = on + f - 1;
        if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end
        s = (sm(t_abs - lag_fr) - sm_mu) / sm_sd;         % z-scored lagged window
        if any(isnan(s)), continue; end
        Sk(end+1,:) = s'; %#ok<AGROW>
        yk(end+1,1) = double(gi==numel(grid) && ~bl.is_censored(b)); %#ok<AGROW>
        tk(end+1,1) = f/fps; %#ok<AGROW>
        lok(end+1,1) = (f-1);     % frames since onset  -> pre-onset  is LEFT  of the onset line %#ok<AGROW>
        ofk(end+1,1) = (dur-f);   % frames until offset -> post-offset is RIGHT of the offset line %#ok<AGROW>
    end
    % pre-entry frames (0..0.5 s), same phase as the grid but DROPPED by left truncation
    grid_tr = allg(allg < entry_fr & allg >= 1);
    Str=[]; ttr=[]; lotr=[]; oftr=[];
    for gi = 1:numel(grid_tr)
        f = grid_tr(gi); t_abs = on + f - 1;
        if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end
        s = (sm(t_abs - lag_fr) - sm_mu) / sm_sd;
        if any(isnan(s)), continue; end
        Str(end+1,:)=s'; ttr(end+1,1)=f/fps; lotr(end+1,1)=(f-1); oftr(end+1,1)=(dur-f); %#ok<AGROW>
    end
    rows_S{k}=Sk; rows_y{k}=yk; rows_tinf{k}=tk; rows_lagonset{k}=lok; rows_offlag{k}=ofk; %#ok<AGROW>
    rows_Str{k}=Str; rows_ttr{k}=ttr; rows_lotr{k}=lotr; rows_oftr{k}=oftr; %#ok<AGROW>
    boutmeta{k} = struct('b',b,'fly',bl.fly(b),'dur',dur,'cens',bl.is_censored(b),'nrow',size(Sk,1)); %#ok<AGROW>
    fprintf('bout %d: fly=%d dur=%d (%.2fs) censored=%d  -> %d at-risk rows\n', ...
        b, bl.fly(b), dur, dur/fps, bl.is_censored(b), size(Sk,1));
end

% ── stack the two bouts with a NaN separator row ─────────────────────────────
sepN = 1;
Sdisp=[]; ydisp=[]; onset_lag=[]; offset_lag=[]; boutRows=zeros(numel(sel),2); truncRows=zeros(numel(sel),2); is_trunc=[]; rlab={};
for k = 1:numel(sel)
    r0 = size(Sdisp,1);
    nTr = size(rows_Str{k},1);  nAr = size(rows_S{k},1);
    Sdisp      = [Sdisp;      rows_Str{k};  rows_S{k}]; %#ok<AGROW>       % truncated rows go ABOVE the at-risk rows
    ydisp      = [ydisp;      zeros(nTr,1); rows_y{k}]; %#ok<AGROW>
    onset_lag  = [onset_lag;  rows_lotr{k}; rows_lagonset{k}]; %#ok<AGROW>
    offset_lag = [offset_lag; rows_oftr{k}; rows_offlag{k}]; %#ok<AGROW>
    is_trunc   = [is_trunc;   true(nTr,1);  false(nAr,1)]; %#ok<AGROW>
    truncRows(k,:) = [r0+1, r0+nTr];
    boutRows(k,:)  = [r0+nTr+1, size(Sdisp,1)];
    for j=1:numel(rows_ttr{k}),  rlab{end+1}=sprintf('%.2f',rows_ttr{k}(j));  end %#ok<AGROW>
    for j=1:numel(rows_tinf{k}), rlab{end+1}=sprintf('%.2f',rows_tinf{k}(j)); end %#ok<AGROW>
    if k < numel(sel)
        Sdisp=[Sdisp; nan(sepN,nLag)]; ydisp=[ydisp; nan(sepN,1)]; onset_lag=[onset_lag; nan(sepN,1)]; offset_lag=[offset_lag; nan(sepN,1)]; is_trunc=[is_trunc; false(sepN,1)]; %#ok<AGROW>
        for j=1:sepN, rlab{end+1}=''; end %#ok<AGROW>
    end
end
nR = size(Sdisp,1);

% ── optional pre-onset masking (within-freeze-only variant) ──────────────────
if mask_preonset
    for r = 1:nR
        if isnan(onset_lag(r)), continue; end
        Sdisp(r, lag_fr > onset_lag(r)) = NaN;   % NaN -> shown as grey (masked to 0)
    end
end

% ════════════════════════════════════════════════════════════════════════
% figure
% ════════════════════════════════════════════════════════════════════════
fh = figure('Color','w','Position',[80 80 650 700]);
tl = tiledlayout(1, 7, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(1, [1, 6])
axM = gca;

% main design-matrix heatmap ------------------------------------------------
hold(axM,'on');
t_rel = -lag_s;                       % past < 0, now = 0, future > 0
[t_ax, ord] = sort(t_rel);            % ascending -> past on the LEFT, now/future on the RIGHT
Sp = Sdisp(:, ord);                   % reorder columns to match the flipped axis
imagesc(axM, t_ax, 1:nR, Sp, 'AlphaData', ~isnan(Sp));
set(axM,'YDir','reverse','Color',[0.86 0.86 0.86]);
cl = max(1, prctile(abs(Sdisp(~isnan(Sdisp))), 99));
colormap(cbrewer2('Reds', []))
caxis(axM, [-3 3]);
xlim(axM,[min(t_ax) max(t_ax)]); ylim(axM,[0.5 nR+0.5]);
xlabel(axM,'lag relative to current frame (s)'); ylabel(axM,'at-risk frame   (time in freeze, s)');

set(axM,'XTick', -kernel_past_s:1:kernel_future_s);
apply_generic(gca, 'font_size', 24)
% thin y labels to ~every 0.5 s: per bout, label the at-risk row nearest each 0.5 s gridline
% (anchor-agnostic — the offset grid need not land on round tinf values)
ytl = repmat({''}, 1, nR);
vv  = str2double(rlab);                       % NaN for '' / separators
for k = 1:numel(sel)
    rr = boutRows(k,1):boutRows(k,2);
    for tg = 0:0.5:ceil(max(vv(rr))*2)/2
        [d,ii] = min(abs(vv(rr) - tg));
        if d <= dt_frames/fps, ytl{rr(ii)} = sprintf('%.2f', vv(rr(ii))); end
    end
end
set(axM,'YTick',1:nR,'YTickLabel',ytl,'TickLength',[0.01 0.01]);   % re-apply AFTER apply_generic (it resets ticks)

% freeze-onset boundary: for each row, lag = (frames since onset)/fps ---------
by = []; bx = [];
for r = 1:nR
    if isnan(onset_lag(r)), by(end+1)=NaN; bx(end+1)=NaN; continue; end %#ok<AGROW>
    by(end+1)=r; bx(end+1)=-onset_lag(r)/fps; %#ok<AGROW>  % t_rel coords: pre-onset is LEFT of this line
end

% freeze-offset boundary: lag to the un-freeze/censor frame = (frames until offset)/fps ------
by2 = []; bx2 = [];
for r = 1:nR
    if isnan(offset_lag(r)), by2(end+1)=NaN; bx2(end+1)=NaN; continue; end %#ok<AGROW>
    by2(end+1)=r; bx2(end+1)=offset_lag(r)/fps; %#ok<AGROW>  % t_rel coords: post-offset is RIGHT of this line
end
plot(axM, bx2, by2, 'w-.', 'LineWidth', 2);
plot(axM, bx, by, 'k-.', 'LineWidth', 2);

% left truncation: fade the pre-entry rows + draw the risk-set entry line ------
xl = xlim(axM);
for k = 1:numel(sel)
    tr = truncRows(k,:);
    if tr(2) < tr(1), continue; end
    patch(axM, [xl(1) xl(2) xl(2) xl(1)], [tr(1)-0.5 tr(1)-0.5 tr(2)+0.5 tr(2)+0.5], ...
        [0.5 0.5 0.5], 'FaceAlpha', 0.9, 'EdgeColor', 'none');
    plot(axM, xl, [tr(2)+0.5 tr(2)+0.5], 'k-', 'LineWidth', 1.5);       % risk-set entry (0.5 s)
    text(axM, xl(1)+0.12, mean(tr), sprintf('Truncated'), ...
        'FontSize', 24, 'Color', 'w', 'VerticalAlignment', 'middle');
end

cb = colorbar(axM); cb.Layout.Tile = 'south';   % attach to the tiledlayout (no manual Position -> no overlap)
cb.Label.String = 'z-scored social motion';
cb.LineWidth = 2;
cb.FontSize = 24;

% outcome strip y ------------------------------------------------------------
nexttile(7)

axY = gca; hold(axY,'on');
imagesc(axY, 1, 1:nR, ydisp, 'AlphaData', ~isnan(ydisp));
set(axY,'YDir','reverse','Color',[0.86 0.86 0.86]);
colormap(axY, [0.93 0.93 0.90; 0.94 0.63 0.15]);  % 0 = pale, 1 = orange event
caxis(axY,[0 1]); xlim(axY,[0.5 1.5]); ylim(axY,[0.5 nR+0.5]);
set(axY,'XTick',1,'XTickLabel',{'y'},'YTick',[],'TickLength',[0 0]);
%title(axY,'Break?','FontSize', 16);
% mark the un-freeze row and the censored bout's final row
for k = 1:numel(sel)
    yb = boutRows(k,:);
    tr = truncRows(k,:);   % fade the truncated rows in the outcome strip too
    if tr(2) >= tr(1)
        patch(axY, [0.5 1.5 1.5 0.5], [tr(1)-0.5 tr(1)-0.5 tr(2)+0.5 tr(2)+0.5], ...
            [0.5 0.5 0.5], 'FaceAlpha', 0.9, 'EdgeColor', 'none');
    end
    if boutmeta{k}.cens
        rectangle(axY,'Position',[0.5 yb(2)-0.5 1 1],'EdgeColor',[0.4 0.4 0.4],'LineStyle','--','LineWidth',1.2);
       % text(axY, 1, yb(2), 'y=0','FontSize', 16,'Color', 'k', 'HorizontalAlignment', 'center');
    else
        %text(axY, 1, yb(2), 'y=1','FontSize', 16, 'Color', 'k', 'HorizontalAlignment', 'center');
    end
end
apply_generic(gca, 'font_size', 24, 'no_y', true, 'no_x', true)

% ── save ─────────────────────────────────────────────────────────────────────
set([axM axY],'Layer','top');   % draw axis box / ticks on top of the image + overlays
axM.Toolbar.Visible='off'; axY.Toolbar.Visible='off';   % keep export clean
if save_out
    outdir = fileparts(mfilename('fullpath'));
    fn = fullfile(outdir, sprintf('design_matrix_example_%s%s.png', grid_anchor, tern(mask_preonset,'_masked','')));
    exportgraphics(fh, fn, 'Resolution', 200);
    [~, base] = fileparts(fn);
    exporter(fh, paths, [base '.pdf']);   % vector PDF -> paths.fig (repo convention)
    fprintf('\nSaved: %s\n       %s\n', fn, fullfile(paths.fig, [base '.pdf']));
end

% ════════════════════════════════════════════════════════════════════════
% local functions
% ════════════════════════════════════════════════════════════════════════
function idx = find_bout(bl, motion_cache, want_cens, dlo, dhi, back_fr, fwd_fr, entry_fr)
% first bout matching censoring flag + duration range with a clean full window
    idx = [];
    for b = 1:height(bl)
        if logical(bl.is_censored(b)) ~= want_cens, continue; end
        dur = bl.dur_frames(b); if dur < dlo || dur > dhi, continue; end
        on = bl.onsets(b); sm = motion_cache(bl.fly(b)); sm = sm(:); L = numel(sm);
        t0 = on + entry_fr - 1;                  % first at-risk frame
        if (t0 - back_fr) < 1 || (on + dur - 1 + fwd_fr) > L, continue; end
        w = sm((t0-back_fr):(on+dur-1+fwd_fr));
        if any(isnan(w)), continue; end
        idx = b; return
    end
end

function C = diverging_bwr(n)
% blue-white-red diverging colormap (no toolbox needed)
    if nargin<1, n=256; end
    m = floor(n/2);
    b = [linspace(0.12,1,m)', linspace(0.35,1,m)', linspace(0.65,1,m)'];   % blue->white
    r = [linspace(1,0.78,n-m)', linspace(1,0.18,n-m)', linspace(1,0.14,n-m)']; % white->red
    C = [b; r];
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
