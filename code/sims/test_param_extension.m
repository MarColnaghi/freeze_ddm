% TEST_PARAM_EXTENSION  Verify the parameter-extended simulation front-end.
%
%   Exercises signal_to_drift.m + simulate_freeze_ddm.m + the x0 extension to
%   sim_leaky_accumulator.m against the checks in the plan:
%     (1) neutral-parameter identity vs the plain drift = sm*beta path
%     (2) exact ReLU/formula port on a hand-computed reference
%     (3) mex vs matlab backend agreement (KS)
%     (4) monotone-response sanity for each knob
%     (5) backward compatibility of the 6-arg sim_leaky_accumulator
%
%   Run headless:
%     matlab -batch "run('sims/test_param_extension.m')"

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'simulators'));
addpath(fullfile(this_dir, 'simulators', 'cpp_mex_code'));

dt    = 1/60;
sigma = 1.0;
theta = 0.35;
beta  = 3.12;
N     = 200000;

% synthetic time-varying, non-negative drive (one trial)
tt = (0:600-1)' * dt;
sm = max(0.5 + 0.4*sin(2*pi*0.2*tt) + 0.1*sin(2*pi*0.7*tt), 0);

npass = 0; nfail = 0;

fprintf('\n=== (2) ReLU/formula port (deterministic) ===\n');
% hand-computed reference of relu_tv_signals_threshold
s_ref  = [0 0.1 0.3 0.6 1.0];
thr = 0.2; off = 0.05; pw = 1.5; eps_p = 1e-6;
relu = max(s_ref - thr, 0);
shifted = (relu > 0) .* (relu + off);
pos = max(shifted, 0);
expected = (pos > 0) .* ((pos + eps_p).^pw - eps_p.^pw);
got = signal_to_drift(s_ref, struct('gain',1,'bias',0, ...
        'relu_threshold',thr,'relu_offset',off,'relu_power',pw));
ok2 = max(abs(got(:) - expected(:))) < 1e-9;
fprintf('  max|got-expected| = %.3e\n', max(abs(got(:) - expected(:))));
print_check('ReLU port matches bayes_fpe formula', ok2);
[npass, nfail] = upd(npass, nfail, ok2);

% neutral preprocessing reduces to gain*signal
d_neutral = signal_to_drift(sm, struct('gain',beta,'bias',0));
ok_neutral = max(abs(d_neutral - sm*beta)) < 1e-12;
print_check('neutral params reduce to drift = beta*sm', ok_neutral);
[npass, nfail] = upd(npass, nfail, ok_neutral);

fprintf('\n=== (1) neutral-parameter identity vs direct mex ===\n');
p0 = struct('dt',dt,'theta',theta,'sigma',sigma,'leak',0,'gain',beta, ...
            'bias',0,'backend','mex','seed',7);
rt_wrap = simulate_freeze_ddm(sm, p0, N);
rt_ref  = sim_ddm_seeded(sm*beta, [dt, sigma^2, 0, theta, 0], N, 7);
a = rt_wrap(~isnan(rt_wrap));  b = rt_ref(~isnan(rt_ref));
ok1 = isequaln(rt_wrap, rt_ref(:));   % same seed + same drift -> identical
fprintf('  median wrap=%.4f  ref=%.4f  cens wrap=%.3f ref=%.3f\n', ...
    median(a), median(b), mean(isnan(rt_wrap)), mean(isnan(rt_ref)));
print_check('wrapper == direct mex call (identical, same seed)', ok1);
[npass, nfail] = upd(npass, nfail, ok1);

fprintf('\n=== (3) mex vs matlab backend agreement (KS) ===\n');
pmex = p0; pmex.backend = 'mex';  pmex.seed = [];
pmat = p0; pmat.backend = 'matlab'; pmat.seed = [];
Nk = 40000;
rng(1); rt_m = simulate_freeze_ddm(sm, pmex, Nk);
rng(1); rt_a = simulate_freeze_ddm(sm, pmat, Nk);
am = rt_m(~isnan(rt_m)); aa = rt_a(~isnan(rt_a));
try
    [~, pval, ks] = kstest2(am, aa);
catch, pval = NaN; ks = NaN; end
ok3 = abs(mean(isnan(rt_m)) - mean(isnan(rt_a))) < 0.01 && ...
      abs(median(am) - median(aa)) < 3*dt;
fprintf('  KS D=%.4f p=%.3g | median mex=%.4f matlab=%.4f | cens %.3f vs %.3f\n', ...
    ks, pval, median(am), median(aa), mean(isnan(rt_m)), mean(isnan(rt_a)));
print_check('mex and matlab backends agree distributionally', ok3);
[npass, nfail] = upd(npass, nfail, ok3);

fprintf('\n=== (4) monotone-response sanity ===\n');
base = struct('dt',dt,'theta',theta,'sigma',sigma,'leak',0,'gain',beta, ...
              'bias',0,'backend','mex','seed',11);
med = @(r) median(r(~isnan(r)));
cens = @(r) mean(isnan(r));

% raising relu_threshold removes evidence -> more censoring
r_lo = simulate_freeze_ddm(sm, setfield(base,'relu_threshold',0.0), N);
r_hi = simulate_freeze_ddm(sm, setfield(base,'relu_threshold',0.5), N);
okA = cens(r_hi) > cens(r_lo);
fprintf('  relu_threshold 0.0->0.5 : censoring %.3f -> %.3f\n', cens(r_lo), cens(r_hi));
print_check('higher relu_threshold increases censoring', okA);
[npass, nfail] = upd(npass, nfail, okA);

% positive delayed_start shifts RTs later by ~dstart
dstart = 0.3;
r_ds = simulate_freeze_ddm(sm, setfield(base,'delayed_start',dstart), N);
okB = (med(r_ds) - med(r_lo)) > 0.5*dstart;
fprintf('  delayed_start=%.2f : median %.4f -> %.4f\n', dstart, med(r_lo), med(r_ds));
print_check('delayed_start shifts RTs later', okB);
[npass, nfail] = upd(npass, nfail, okB);

% positive sensory_delay (delays evidence) -> later / more censored
r_sd = simulate_freeze_ddm(sm, setfield(base,'sensory_delay',0.3), N);
okC = med(r_sd) >= med(r_lo) - dt || cens(r_sd) > cens(r_lo);
fprintf('  sensory_delay=0.3 : median %.4f (base %.4f) cens %.3f (base %.3f)\n', ...
    med(r_sd), med(r_lo), cens(r_sd), cens(r_lo));
print_check('positive sensory_delay does not speed up responses', okC);
[npass, nfail] = upd(npass, nfail, okC);

% start point x0 = frac*theta > 0 speeds up crossings
r_x0 = simulate_freeze_ddm(sm, setfield(base,'initial_value_frac',0.5), N);
okD = med(r_x0) < med(r_lo);
fprintf('  initial_value_frac=0.5 : median %.4f -> %.4f\n', med(r_lo), med(r_x0));
print_check('positive start point speeds up crossings', okD);
[npass, nfail] = upd(npass, nfail, okD);

% contaminant lapses replace ~pcont of RTs with Uniform(0,t_max): here normal
% crossings are ~0.15s while t_max~10s, so lapses inject late-RT mass. Expect
% frac(rt>1s) ~ pcont*(1 - 1/t_max).
pcont = 0.2;
r_ct = simulate_freeze_ddm(sm, setfield(base,'contaminant_prob',pcont), N);
t_max = numel(sm) * dt;
frac_late_base = mean(r_lo(~isnan(r_lo)) > 1.0);
frac_late_ct   = mean(r_ct(~isnan(r_ct)) > 1.0);
expected_late  = pcont * (1 - 1.0 / t_max);
okE = frac_late_base < 0.02 && abs(frac_late_ct - expected_late) < 0.03;
fprintf('  contaminant_prob=%.1f : frac(rt>1s) %.3f -> %.3f (expected ~%.3f)\n', ...
    pcont, frac_late_base, frac_late_ct, expected_late);
print_check('contaminant lapses inject ~pcont uniform-RT mass', okE);
[npass, nfail] = upd(npass, nfail, okE);

fprintf('\n=== (4b) start point: initial_value_frac vs x0 ===\n');
% the two conventions must agree when they denote the same start point
frac = 0.4;
p_frac = base; p_frac.initial_value_frac = frac;
p_abs  = base; p_abs.x0 = frac * theta;
r_frac = simulate_freeze_ddm(sm, p_frac, N);
r_abs  = simulate_freeze_ddm(sm, p_abs,  N);
okF = isequaln(r_frac, r_abs);   % same seed + same resolved x0 -> identical
fprintf('  frac=%.2f (x0=%.3f): median frac=%.4f abs=%.4f\n', ...
    frac, frac*theta, med(r_frac), med(r_abs));
print_check('initial_value_frac and equivalent x0 give identical RTs', okF);
[npass, nfail] = upd(npass, nfail, okF);

% x0 alone resolves correctly
[~, info_abs] = simulate_freeze_ddm(sm, p_abs, 1);
okG = abs(info_abs.x0 - frac*theta) < 1e-12;
print_check('p.x0 is used verbatim (absolute units)', okG);
[npass, nfail] = upd(npass, nfail, okG);

% giving both must error
p_both = base; p_both.initial_value_frac = 0.4; p_both.x0 = 1.0;
okH = false;
try
    simulate_freeze_ddm(sm, p_both, 10);
catch ME
    okH = strcmp(ME.identifier, 'simulate_freeze_ddm:AmbiguousStartPoint');
    fprintf('  errored as expected: %s\n', ME.identifier);
end
print_check('specifying both initial_value_frac and x0 errors', okH);
[npass, nfail] = upd(npass, nfail, okH);

% neither -> x0 = 0
[~, info_none] = simulate_freeze_ddm(sm, base, 1);
okI = info_none.x0 == 0;
print_check('neither given -> x0 = 0', okI);
[npass, nfail] = upd(npass, nfail, okI);

fprintf('\n=== (6) delay equivalence vs bayes_fpe + pre-onset stream ===\n');
% (6a) THE reference check: our drift must equal bayes_fpe's solver input, i.e.
% delay_tv_signals_with_zero_pad(raw, eff=sensory_delay-delayed_start) followed by
% the delayed_start slice. Reimplemented here from tv_core.py:326 and :383-400.
Tp = 200; n_pre = 48;                       % 48 frames = 0.8 s of pre-onset margin
rng(3); full_raw = max(cumsum(randn(Tp+n_pre,1))*0.05 + 1, 0);
post_raw = full_raw(n_pre+1:end);
st_post  = (0:Tp-1)'*dt;
ref_zero_pad = @(sig, lag, st) ref_delay(sig, lag, st);

worst = 0;
for Sd = [0 0.15 0.3 0.5]
  for Dd = [0 0.15 0.3 0.5]
    eff  = Sd - Dd;                          % effective_solver_sensory_delay
    ref  = ref_zero_pad(post_raw, eff, st_post);
    d    = signal_to_drift(post_raw, struct('dt',dt,'sensory_delay',Sd,'gain',1,'bias',0));
    mine = d(round(Dd/dt)+1:end);            % delayed_start slice
    worst = max(worst, max(abs(mine - ref(1:numel(mine)))));
  end
end
ok6a = worst < 1e-9;
fprintf('  worst |drift - bayes_fpe| over 4x4 (sensory_delay x delayed_start) = %.3e\n', worst);
print_check('drift matches bayes_fpe delay pipeline (n_pre=0)', ok6a);
[npass, nfail] = upd(npass, nfail, ok6a);

% (6b) zero pre-onset must reduce EXACTLY to the zero-pad case
L = 0.4;
d_pad  = signal_to_drift(post_raw, struct('dt',dt,'sensory_delay',L,'gain',1,'bias',0));
d_zero = signal_to_drift([zeros(n_pre,1); post_raw], ...
             struct('dt',dt,'sensory_delay',L,'gain',1,'bias',0,'n_pre',n_pre));
ok6b = isequal(size(d_pad), size(d_zero)) && max(abs(d_pad - d_zero)) < 1e-12;
print_check('zero pre-onset reduces exactly to the zero-pad case', ok6b);
[npass, nfail] = upd(npass, nfail, ok6b);

% (6c) real pre-onset IS consumed: drift on [0,L) equals raw[tau-L], not zero
d_real = signal_to_drift(full_raw, ...
             struct('dt',dt,'sensory_delay',L,'gain',1,'bias',0,'n_pre',n_pre));
nL     = round(L/dt);
expect = full_raw(n_pre+1-nL : n_pre-nL+nL);
ok6c   = max(abs(d_real(1:nL) - expect)) < 1e-9 && max(abs(d_pad(1:nL))) < 1e-12;
fprintf('  on [0,L): real mean=%.4f vs zero-pad max=%.3e\n', ...
    mean(d_real(1:nL)), max(abs(d_pad(1:nL))));
print_check('real pre-onset stream is used on [0,sensory_delay)', ok6c);
[npass, nfail] = upd(npass, nfail, ok6c);

% (6d) beyond the supplied margin, fall back to zero-pad
Lbig  = (n_pre+20)*dt;
d_big = signal_to_drift(full_raw, ...
             struct('dt',dt,'sensory_delay',Lbig,'gain',1,'bias',0,'n_pre',n_pre));
ok6d  = all(abs(d_big(1:20)) < 1e-12) && any(abs(d_big(25:60)) > 0);
print_check('falls back to zero-pad beyond the supplied margin', ok6d);
[npass, nfail] = upd(npass, nfail, ok6d);

% (6e) trimming contract: output always aligned to t=0 and the same length
ok6e = numel(d_pad) == Tp && numel(d_real) == Tp;
[~, fl] = signal_to_drift(full_raw, ...
              struct('dt',dt,'sensory_delay',L,'gain',1,'bias',0,'n_pre',n_pre));
ok6e = ok6e && abs(fl.t(n_pre+1)) < 1e-12 && numel(fl.t) == Tp + n_pre;
print_check('n_pre trims to onset; full.t has t=0 at n_pre+1', ok6e);
[npass, nfail] = upd(npass, nfail, ok6e);

% (6f) scientific sanity: real pre-onset slows RTs LESS than zero-pad
pz = base; pz.sensory_delay = L;
pr = base; pr.sensory_delay = L; pr.n_pre = n_pre;
sm_pre = [sm(1)*ones(n_pre,1); sm];          % plausible non-zero pre-onset stream
rt_z = simulate_freeze_ddm(sm,     pz, N);
rt_r = simulate_freeze_ddm(sm_pre, pr, N);
ok6f = med(rt_r) < med(rt_z);
fprintf('  sensory_delay=%.2f: median zero-pad=%.4f vs real pre-onset=%.4f\n', ...
    L, med(rt_z), med(rt_r));
print_check('real pre-onset slows RTs less than zero-pad', ok6f);
[npass, nfail] = upd(npass, nfail, ok6f);

fprintf('\n=== (5) backward compatibility of sim_leaky_accumulator ===\n');
rng(3); k6 = sim_leaky_accumulator(sm*beta, theta, sigma, 0, dt);          % 5-arg
rng(3); k7 = sim_leaky_accumulator(sm*beta, theta, sigma, 0, dt, [], 0);   % 7-arg x0=0
ok5 = isequaln(k6, k7);
print_check('6/7-arg calls with x0=0 unchanged', ok5);
[npass, nfail] = upd(npass, nfail, ok5);

fprintf('\n==================== %d passed, %d failed ====================\n', npass, nfail);
if nfail > 0, error('test_param_extension: %d checks failed', nfail); end

% ---- helpers ----
function print_check(name, cond)
    if cond, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  [%s] %s\n', tag, name);
end
function [np, nf] = upd(np, nf, cond)
    if cond, np = np + 1; else, nf = nf + 1; end
end
function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
function out = ref_delay(sig, lag, st)
% Reference reimplementation of bayes_fpe delay_tv_signals_with_zero_pad
% (tv_core.py:326): sample at t-lag, left edge -> 0, right edge -> last value.
    tq  = st - lag;
    out = interp1(st, sig, tq, 'linear', NaN);
    out(tq < st(1))   = 0;
    out(tq > st(end)) = sig(end);
    out = out(:);
end
