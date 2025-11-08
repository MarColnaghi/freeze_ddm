function [nll, f, fd, F] = nll_fly_ddm_newer(params, bouts, points, model_num, iid, plot_flag, extra)

% fd = freeze_durations
% y1 = average_motion post to the freeze onset
% y2 = average_motion previous to the freeze onset
% y3 = speed_focal
% y4 = loom_number
% y5 = speed of the loom

if nargin < 7
    extra = [];
end

model_func = str2func(model_num);
pattern = '_(\D+)\d+';
tok = regexp(model_num, pattern, 'tokens');

if strcmp(plot_flag, 'p')

    tbl = table();
    tbl.durations_s = [0:1/300:(points.censoring + 2)]';
    fd = tbl.durations_s;
    f = zeros(height(tbl.durations_s), 1);
    F = zeros(height(tbl.durations_s), 1);
    G = zeros(height(tbl.durations_s), 1);

    for idx_bout = 1:height(bouts)
        tbl.sm = bouts.sm(idx_bout)*ones(height(tbl),1);
        tbl.smp = bouts.smp(idx_bout)*ones(height(tbl),1);
        tbl.fs = bouts.fs(idx_bout)*ones(height(tbl),1);
        tbl.ln = bouts.ln(idx_bout)*ones(height(tbl),1);
        tbl.ls = bouts.ls(idx_bout)*ones(height(tbl),1);
        tbl.intercept = bouts.intercept(idx_bout)*ones(height(tbl),1);

        % ec.soc_mot_array = extra.soc_mot_array(idx_bout, :);
        g = comp_loglikelihood(params, tbl, points, model_func, iid, tok, extra);
        %[~, G] = comp_loglikelihood(params, tbl, points, model_func, iid, tok, extra);

        F = F + G;
        f = f + exp(g);
        trapz(tbl.durations_s(tbl.durations_s > points.truncation & tbl.durations_s <= points.censoring), f(tbl.durations_s > points.truncation &  tbl.durations_s <= points.censoring));
    end

    f = f ./ height(bouts);
    F = F ./ height(bouts);
    trapz(tbl.durations_s(tbl.durations_s > points.truncation & tbl.durations_s <= points.censoring), f(tbl.durations_s > points.truncation & tbl.durations_s <= points.censoring)) + f(end)

end
arr = unique(bouts.fly)';
if strcmp(iid, 'dep')
    g = zeros(1, max(arr));

    for idx_flies = arr
        idx = bouts.fly == idx_flies;
        g(idx_flies) = sum(comp_loglikelihood(params, bouts(idx, :), points, model_func, iid, tok, extra));
    end

elseif strcmp(iid, 'iid')
    g = zeros(1, height(bouts));
    g = comp_loglikelihood(params, bouts, points, model_func, iid, tok, extra);

end

nll = -sum(g);

end

function [log_g] = comp_loglikelihood(x, bouts_individual_fly, points, model_func, iid, tok, extra)

bif = bouts_individual_fly;

ts = bif.durations_s;
y = table;
y.sm = bif.sm;
y.smp = bif.smp;
y.fs = bif.fs;
y.ln = bif.ln;
y.ls = bif.ls;
y.intercept = bif.intercept;

model = model_func();
[gt, lbl] = get_ground_truth_vector(model);
lbl = lbl(~isnan(gt));
gt_table = array2table(x, 'VariableNames', lbl);
out = evaluate_model(model, gt_table, y);

if ~isfield(model, 'tndt')
    out.tndt = zeros(size(out, 1), 1);
end

% assess freeze duration categories
bet = bif.durations_s <= points.censoring;
abo = bif.durations_s > points.censoring;

g = zeros(size(ts));

if strcmp('iid', iid)

    if  strcmp('sddm', tok{1})

        bet = bif.durations_s > out.tndt & bif.durations_s <= points.censoring;

        % pdf and cdf for single bound ddm
        [pdf, cdf] = pdf_cdf({'ddm'});

        f = @(ts, inds) pdf.ddm(ts, out.mu(inds), out.theta(inds), out.tndt(inds));
        F = @(ts, inds) cdf.ddm(ts, out.mu(inds), out.theta(inds), out.tndt(inds));

        if ~isempty(points.truncation) && points.truncation > out.tndt(1)
            trunc_factor = @(inds) 1 - F(points.truncation, inds);
        else
            trunc_factor = @(inds) ones(size(ts(inds)));
        end

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = (1 - F(points.censoring, abo)) ./ trunc_factor(abo);

        g = max(g, 1e-5);
        log_g = log(g);

    elseif  strcmp('dddm', tok{1})

        bet = bif.durations_s >= out.tndt & bif.durations_s <= points.censoring;
        abo = bif.durations_s > points.censoring;

        % pdf and cdf for single bound ddm
        [pdf, cdf] = pdf_cdf({'ddm','kde'});

        pdf_ddm_raw = pdf.ddm;
        cdf_ddm_raw = cdf.ddm;

        pdf.ddm = @(ts, mu, theta, ndt) guard_ddm(pdf_ddm_raw, ts, mu, theta, ndt);
        cdf.ddm = @(ts, mu, theta, ndt) guard_ddm(cdf_ddm_raw, ts, mu, theta, ndt);

        f = @(ts, inds) out.pmix(inds) .* pdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds)) + ...
            (1 - out.pmix(inds)) .* pdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));
        F = @(ts, inds) out.pmix(inds) .* cdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds)) + ...
            (1 - out.pmix(inds)) .* cdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));

        t0   = points.truncation;
        C    = points.censoring;
        epsN = 1e-12;

        % One consistent truncation factor: 1 - F_mix(t0) per index
        trunc_factor = @(inds) max(1 - F(t0, inds), epsN);

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = (1 - F(points.censoring, abo)) ./ trunc_factor(abo);

        g = max(g, 1e-5);
        log_g = log(g);

    elseif  strcmp('ksddm', tok{1})

        below = bif.durations_s <  out.tndt;
        bet   = bif.durations_s >= out.tndt & bif.durations_s <= points.censoring;
        abo   = bif.durations_s >  points.censoring;

        % pdf and cdf (UNtruncated) for the two components
        [pdf, cdf] = pdf_cdf({'ddm','kde'});

        pdf_ddm_raw = pdf.ddm;
        cdf_ddm_raw = cdf.ddm;

        pdf.ddm = @(ts, mu, theta, ndt) guard_ddm(pdf_ddm_raw, ts, mu, theta, ndt);
        cdf.ddm = @(ts, mu, theta, ndt) guard_ddm(cdf_ddm_raw, ts, mu, theta, ndt);

        % Mixers
        f_kde = @(ts, inds, extra) (1 - out.pmix(inds)) .* pdf.kde(ts, extra);
        F_kde = @(t,  inds, extra) (1 - out.pmix(inds)) .* cdf.kde(t,  extra);

        f_ddm = @(ts, inds) out.pmix(inds) .* pdf.ddm(ts, out.mu(inds), out.theta(inds), out.tndt(inds));
        F_ddm = @(t,  inds) out.pmix(inds) .* cdf.ddm(t,  out.mu(inds), out.theta(inds), out.tndt(inds));

        f = @(ts, inds, extra) f_ddm(ts, inds) + f_kde(ts, inds, extra);
        F = @(t,  inds, extra) F_ddm(t,  inds) + F_kde(t,  inds, extra);

        t0   = points.truncation;
        C    = points.censoring;
        epsN = 1e-12;

        % One consistent truncation factor: 1 - F_mix(t0) per index
        trunc_factor = @(inds) max(1 - F(t0, inds, extra), epsN);

        % Likelihoods
        g          = nan(size(bif.durations_s));
        g(below)   = f_kde(ts(below), below, extra) ./ trunc_factor(below);
        g(bet)     = f(ts(bet),   bet,   extra)     ./ trunc_factor(bet);
        g(abo)     = (1 - F(C, abo, extra))         ./ trunc_factor(abo);

        g      = max(g, 1e-5);
        log_g  = log(g);

    elseif  strcmp('exp', tok{1})

        % pdf and cdf for single bound ddm
        [pdf, cdf] = pdf_cdf({'ddm', 'exp'});

        f = @(ts, inds) pdf.exp(ts, out.lambda(inds));
        F = @(ts, inds) cdf.exp(ts, out.lambda(inds));

        if ~isempty(points.truncation)
            trunc_factor = @(inds) 1 - F(points.truncation, inds);
        else
            trunc_factor = @(inds) ones(size(ts(inds)));
        end

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = (1 - F(points.censoring, abo)) ./ trunc_factor(abo);

        g      = max(g, 1e-5);
        log_g  = log(g);

    elseif  strcmp('ed', tok{1})

        fs = 60;
        % out.tndt = ones(length(out.theta), 1) * 0.15;

        below = bif.durations_s <  out.tndt;
        bet   = bif.durations_s >= out.tndt & bif.durations_s <= points.censoring;
        abo   = bif.durations_s >  points.censoring;

        if size(extra.soc_mot_array, 1) == 1
            out.mu = repmat(extra.soc_mot_array, height(out.theta), 1) .* x(1) .* (1/fs);
        else
            out.mu = extra.soc_mot_array .* x(1) .* (1/fs);
        end

        [pdf, cdf] = pdf_cdf({'ed'});

        f = @(ts, inds) pdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);
        F = @(ts, inds) cdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);

        if ~isempty(points.truncation)
            trunc_factor = @(inds) 1 - F(points.truncation, inds);
        else
            trunc_factor = @(inds) ones(size(ts(inds)))';
        end

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = F(points.censoring, abo) ./ trunc_factor(abo);

        g      = max(g, 1e-5);
        log_g  = log(g);

%         log_g = g;
%         ndt_prior: [1 x n_ndt] probability distribution over ndt values

%         fs = 60;
% 
%         if size(extra.soc_mot_array, 1) == 1
%             out.mu = repmat(extra.soc_mot_array, height(out.theta), 1) .* x(1) .* (1/fs);
%         else
%             out.mu = extra.soc_mot_array .* x(1) .* (1/fs);
%         end
% 
%         n_ndt = 31;
%         ndt_values = (0:(n_ndt-1)) / fs; % convert frames to seconds
%         ndt_prior = ones(1, n_ndt) / n_ndt;
% 
%         [pdf, cdf] = pdf_cdf({'ed'});
% 
%         n_trials = length(ts);
%         log_liks = zeros(n_trials, n_ndt);
% 
%         for i = 1:n_ndt
% 
%             current_ndt = ndt_values(i);
% 
%             out.tndt = current_ndt *  ones(length(out.theta), 1);
%             
%             bet = ts > out.tndt & ts - out.tndt < points.censoring;
%             abo = ts - out.tndt >= points.censoring;
% 
%             f = @(ts, inds) pdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);
%             F = @(ts, inds) cdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);
% 
%             g(bet) = f(ts(bet), bet);
%             g(abo) = F(points.censoring, abo);
%             pdf_vals = g;
%             log_liks(:, i) = log(pdf_vals(:) + 1e-10);
%         end
% 
%         % Marginal: sum_ndt p(data|params,ndt) * p(ndt)
%         % In log space with logsumexp trick
% 
%         max_ll = max(log_liks, [], 2);  % [n_trials × 1]
%         trial_marginal_ll = max_ll + log(sum(exp(log_liks - max_ll) .* ndt_prior, 2));
% 
%         log_g = trial_marginal_ll;
    end
end

end

function y = guard_ddm(fun, t, mu, th, ndt)
    % Returns 0 for entries where t <= ndt; calls `fun` otherwise
    mu  = mu(:); th = th(:); ndt = ndt(:);
    if isscalar(t)
        t = repmat(t, size(mu));     % handle scalar t0/C with vector params
    else
        t = t(:);
    end
    y = zeros(size(mu));
    m = t > ndt;
    if any(m)
        y(m) = fun(t(m), mu(m), th(m), ndt(m));
    end
end
