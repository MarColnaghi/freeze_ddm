function nullres = sm_circshift_null(bl, motion_cache, opts, n_shuffle, min_shift_s)
% SM_CIRCSHIFT_NULL  Circular-time-shift null for the social-motion hazard kernel.
%
%   nullres = sm_circshift_null(bl, motion_cache, opts, n_shuffle, min_shift_s)
%
%   Tests whether the recovered hazard kernel reflects genuine sm -> escape
%   coupling rather than just sm autocorrelation + the timing of escapes. Each
%   shuffle circularly rotates EVERY fly's sm trace by an independent random
%   offset before re-reading the kernel windows on the SAME bouts. The rotation
%   preserves each trace's full autocorrelation / amplitude distribution and
%   breaks only its alignment to the real un-freeze times, so the shuffled
%   kernels are the null expected from autocorrelation structure alone.
%
%   This is a thin wrapper: run_hazard_kernel does the rotation INTERNALLY via
%   shuffle_mode='circshift' (design built once, each shuffle only re-indexes the
%   shifted windows + a warm-started refit at the real-fit lambda; the leak-fit /
%   2000x bootstrap is skipped). The wrap seam is blanked to the mean for any
%   window that straddles it (see run_hazard_kernel; reported as seam_frac).
%
%   n_shuffle (default 200); min_shift_s (default 5 s) the minimum rotation.
%
%   nullres fields:
%       t_rel, caus_mask              lag axis (from the real fit)
%       kernel_real                   real recovered kernel
%       kernels_null  (nLag x n_shuffle)
%       null_lo/null_hi/null_med      2.5 / 97.5 / 50 pct band over shuffles, per lag
%       peak_real, sum_real           causal peak / sum of the real kernel
%       peak_null, sum_null  (n_shuffle x 1)
%       p_peak, p_sum                 one-sided (frac null >= real, +1 smoothed)
%       n_events, seam_frac, n_shuffle, lambda

    if nargin < 4 || isempty(n_shuffle),   n_shuffle   = 200; end
    if nargin < 5 || isempty(min_shift_s), min_shift_s = 5;   end

    o = opts;
    o.shuffle_mode     = 'circshift';
    o.n_shuffle        = n_shuffle;
    o.circ_min_shift_s = min_shift_s;
    o.do_leakfit = false;     % null only needs the kernel
    o.do_unpen   = false;
    o.do_stage2  = false;
    o.do_plot    = false;

    res = run_hazard_kernel(bl, motion_cache, o);

    t_rel = res.t_rel; cm = res.caus_mask;
    kernel_real  = res.kernel;
    kernels_null = res.kernels_sh;                 % nLag x n_shuffle

    peak_real = max(kernel_real(cm));  sum_real = sum(kernel_real(cm));
    peak_null = max(kernels_null(cm,:), [], 1).';
    sum_null  = sum(kernels_null(cm,:), 1).';
    p_peak = (1 + sum(peak_null >= peak_real)) / (1 + n_shuffle);
    p_sum  = (1 + sum(sum_null  >= sum_real )) / (1 + n_shuffle);

    null_lo  = prctile(kernels_null, 2.5, 2);
    null_hi  = prctile(kernels_null, 97.5,2);
    null_med = prctile(kernels_null, 50,  2);

    nullres = struct('t_rel',t_rel,'caus_mask',cm, ...
        'kernel_real',kernel_real,'kernels_null',kernels_null, ...
        'null_lo',null_lo,'null_hi',null_hi,'null_med',null_med, ...
        'peak_real',peak_real,'sum_real',sum_real, ...
        'peak_null',peak_null,'sum_null',sum_null,'p_peak',p_peak,'p_sum',p_sum, ...
        'n_events',res.n_events,'seam_frac',getf(res,'seam_frac',NaN), ...
        'n_shuffle',n_shuffle,'lambda',res.lambda);

    fprintf('circshift null | events=%d | seam-blanked %.1f%% of windows | peak p=%.3f, sum p=%.3f\n', ...
        res.n_events, 100*nullres.seam_frac, p_peak, p_sum);

    % ── plot: real kernel vs null band ─────────────────────────────────────────
    if geto(opts,'do_plot',true)
        figure('Color','w','Position',[100 100 760 470]); hold on
        fill([t_rel;flipud(t_rel)],[null_hi;flipud(null_lo)],[.6 .6 .6], ...
            'FaceAlpha',.3,'EdgeColor','none','DisplayName','null 95% (circshift)')
        plot(t_rel,null_med,    'Color',[.4 .4 .4],'LineWidth',1,'DisplayName','null median')
        plot(t_rel,kernel_real, 'Color',[.2 .5 .8],'LineWidth',2,'DisplayName','real kernel')
        yl = ylim; tac = max(t_rel);
        hp = patch([0 tac tac 0],[yl(1) yl(1) yl(2) yl(2)],[.92 .92 .92], ...
            'EdgeColor','none','HandleVisibility','off'); uistack(hp,'bottom'); ylim(yl);
        yline(0,'k:'); xline(0,'--k','escape');
        xlabel('Time relative to un-freeze (s)   (causal \tau\leq0 | acausal \tau>0 shaded)')
        ylabel('Kernel \beta(\tau)')
        title(sprintf('Hazard kernel vs circular-shift null   (peak p=%.3f, \\Sigma p=%.3f)', p_peak, p_sum))
        legend('Location','northwest','box','off')
    end
end

% ════════════════════════════════════════════════════════════════════════
function v = geto(s,f,d), if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end, end
function v = getf(s,f,d), if isfield(s,f), v=s.(f); else, v=d; end, end
