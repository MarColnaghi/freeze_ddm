% ---------- settings ----------
fps = 60;                          % social-motion sampling rate
dt  = 0.1;                         % analysis bin (s) -> averages 6 frames/bin
K   = 30;                          % lags -> K*dt = 3.0 s of history
floor_s = 0.5;                     % detection floor: shorter freezes are absent (left-truncation)
k0 = round(floor_s/dt) + 1;        % first detectable bin (delayed entry)
flyidx = findgroups(bouts_proc.fly);

% ---------- motion traces, downsampled 60 Hz -> dt bins ----------
motion_ts = extract_sm_from_bouts(bouts_proc, 'type','onlyfreeze', ...
    'output_type','cell', 'norm_factor',10, 'cache','motion_cache');
factor = max(1, round(dt*fps));    % frames per bin (= 6 for dt = 0.1)
ds = @(v) mean(reshape(v(1:floor(numel(v)/factor)*factor), factor, []), 1).';

% ---------- build the bout-period table (at-risk only from the 0.5 s floor) ----------
buf = cell(height(bouts_proc),1);
for i = 1:height(bouts_proc)
    m  = ds(motion_ts{i}(:));                               % motion at dt resolution
    ni = min(numel(m), round(bouts_proc.durations_s(i)/dt));
    ev = double(bouts_proc.contacts(i) == 0);               % 1 = resumption, 0 = collision
    last = ni - (1 - ev);                                   % censored: survives to ni-1
    if last < k0, continue; end                             % must reach the detectable floor
    k = (k0:last).';                                        % delayed entry: rows from 0.5 s on
    L = zeros(numel(k), K+1);
    for j = 0:K
        idx = k - j; ok = idx >= 1;                         % lags still use [0,0.5) real motion
        L(ok, j+1) = m(idx(ok));
    end
    buf{i} = [(k-1)*dt, double(ev & k==ni), repmat(double(flyidx(i)),numel(k),1), L];
end
buf = buf(~cellfun(@isempty, buf));
P = array2table(cell2mat(buf), 'VariableNames', cellstr(["t","y","fly","m"+(0:K)]));

% ---------- fit: full (history) vs current-only ----------
lagterms = strjoin("m"+(0:K), " + ");
full = fitglm(P, "y ~ t + t^2 + " + lagterms, 'Distribution','binomial');
now_ = fitglm(P, "y ~ t + t^2 + m0",          'Distribution','binomial');
LR = now_.Deviance - full.Deviance;
fprintf('history LR: chi2 = %.1f, df = %d, p = %.2e\n', LR, K, 1 - chi2cdf(LR, K));

% ---------- kernel ----------
figure
kernel = full.Coefficients.Estimate(startsWith(full.CoefficientNames,'m'));
plot((0:K)*dt, kernel, '-o', 'LineWidth', 2)
xlabel('lag into the past (s)'); ylabel('kernel weight (log-odds / unit motion)')

%%
dt = 0.1; K = 15;                         % trim window to 1.5 s; the structure is all < 1 s
nb = 6;                                   % raised-cosine bumps
lags = (0:K)'; ctrs = linspace(0,K,nb); w = K/(nb-1);
Bspl = 0.5*(1+cos(max(min((lags-ctrs)/w*pi,pi),-pi)));   % (K+1) x nb basis

buf = cell(height(bouts_proc),1);
for i = 1:height(bouts_proc)
    m  = ds(motion_ts{i}(:));
    ni = min(numel(m), round(bouts_proc.durations_s(i)/dt));
    ev = double(bouts_proc.contacts(i) == 0);
    last = ni - (1 - ev); if last < 1, continue; end
    k = (1:last)'; L = zeros(last, K+1);
    for j = 0:K
        idx = k - j; ok = idx >= 1; L(ok, j+1) = m(idx(ok));
    end
    buf{i} = [(k-1)*dt, double(ev & k==ni), repmat(double(flyidx(i)),last,1), L*Bspl];
end
buf = buf(~cellfun(@isempty,buf));
P = array2table(cell2mat(buf), 'VariableNames', cellstr(["t","y","fly","x"+(1:nb)]));

xterms = strjoin("x"+(1:nb), " + ");
full = fitglm(P, "y ~ t + t^2 + " + xterms, 'Distribution','binomial');

bsel  = startsWith(full.CoefficientNames,'x');
beta  = full.Coefficients.Estimate(bsel);
Cov   = full.CoefficientCovariance(bsel,bsel);
kernel = Bspl*beta;  se = sqrt(diag(Bspl*Cov*Bspl'));     % kernel + analytic CI
errorbar((0:K)*dt, kernel, 1.96*se, 'LineWidth', 2)
xlabel('lag into the past (s)'); ylabel('kernel weight (log-odds / unit motion)')

%%

% 2) how collinear are the regressors / how fast does motion decorrelate?
allm = cell2mat(cellfun(@(v) ds(v(:)), motion_ts, 'UniformOutput', false));
ac = xcorr(allm-mean(allm), K, 'coeff'); ac = ac(K+1:end);
fprintf('motion autocorr at 0.1 / 0.3 s: %.2f / %.2f   cond(X) = %.0f\n', ...
        ac(2), ac(4), cond(P{:, "x"+(1:nb)}));