function [nll, f, fd] = nll_fly_ddm_newer(params, bouts, points, model_num, iid, plot_flag, extra)

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

    
    fd = [1/60:1/60:(points.censoring + 1/60)]';
    num_times = height(fd);
    num_bouts = height(bouts);

    tbl = table();
    tbl.durations_s = fd;
    tbl.sm = zeros(num_times, 1);
    tbl.smp = zeros(num_times, 1);
    tbl.fs = zeros(num_times, 1);
    tbl.ln = zeros(num_times, 1);
    tbl.ls = zeros(num_times, 1);
    tbl.intercept = zeros(num_times, 1);
    
    f = zeros(num_times, 1);
    
    hold on

    for idx_bout = 1:num_bouts

        tbl.sm(:) = bouts.sm(idx_bout);
        tbl.smp(:) = bouts.smp(idx_bout);
        tbl.fs(:) = bouts.fs(idx_bout);
        tbl.ln(:) = bouts.ln(idx_bout);
        tbl.ls(:) = bouts.ls(idx_bout);
        tbl.intercept(:) = bouts.intercept(idx_bout);

        if strcmp('ed', tok{1}) || strcmp('ded', tok{1})
            ec.soc_mot_array = extra.soc_mot_array(idx_bout, :);
        else
            ec = [];
        end
        g = comp_loglikelihood(params, tbl, points, model_func, iid, tok, ec);

        if ~isempty(points.truncation)
            % fprintf('bouts %d: sum(f): = %d \n', idx_bout, trapz(fd(fd >= points.truncation & fd <= points.censoring) , exp(g(fd >= points.truncation & fd <= points.censoring))) + exp(g(end)))
            fprintf('bouts %d: sum(f): = %d \n', idx_bout, sum(exp(g(fd > points.truncation & fd <= points.censoring))) + exp(g(end)))

        else
            fprintf('bouts %d: sum(f): = %d \n', idx_bout, trapz(exp(g(fd > 0 & fd <= points.censoring))) + exp(g(end)))
        end

        f = f + exp(g);

    end

    f = f ./ num_bouts;
    nll = [];
else

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
if isfield(extra, 'tndt')
    model = rmfield(model, 'tndt');
end

[gt, lbl] = get_ground_truth_vector(model);
lbl = lbl(~isnan(gt));
gt_table = array2table(x, 'VariableNames', lbl);

if strcmp('ed', tok{1}) || strcmp('ded', tok{1})
    if isfield(extra, 'soc_mot_array')
        if size(extra.soc_mot_array, 1) == 1
            y.sm = repmat(extra.soc_mot_array, height(y), 1);
        else
            y.sm = extra.soc_mot_array;
        end

    end
end

out = evaluate_model(model, gt_table, y);

if ~isfield(model, 'tndt')
    out.tndt = zeros(size(out, 1), 1);
    if isfield(extra, 'tndt')
        out.tndt = ones(size(out, 1), 1) * extra.tndt;
    end
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

        epsN = 1e-12;
        g      = max(g, epsN);
        log_g  = log(g);

    elseif  strcmp('ed', tok{1})

        fs = 60;

        below = bif.durations_s <=  out.tndt;
        bet   = bif.durations_s > out.tndt & bif.durations_s <= points.censoring;
        abo   = bif.durations_s >  points.censoring;

        [pdf, cdf] = pdf_cdf({'ed'});

        f = @(ts, inds) pdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);
        F = @(ts, inds) cdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);

        epsN = 1e-12;

        if ~isempty(points.truncation)
            trunc_factor = @(inds) max(F(points.truncation, inds), epsN) ;
        else
            trunc_factor = @(inds) ones(size(ts(inds)))';
        end

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = F(points.censoring, abo) ./ trunc_factor(abo);

        g      = max(g, epsN);
        log_g  = log(g);

    elseif  strcmp('simed', tok{1})

        fs = 60;

        out.mu = zeros(height(out), 1) .* ones(height(out), 630);

        below = bif.durations_s <=  out.tndt;
        bet   = bif.durations_s > out.tndt & bif.durations_s <= points.censoring;
        abo   = bif.durations_s >  points.censoring;

        [pdf, cdf] = pdf_cdf({'ed'});

        f = @(ts, inds) pdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);
        F = @(ts, inds) cdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);

        epsN = 1e-12;

        if ~isempty(points.truncation)
            trunc_factor = @(inds) max(F(points.truncation, inds), epsN) ;
        else
            trunc_factor = @(inds) ones(size(ts(inds)))';
        end

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = F(points.censoring, abo) ./ trunc_factor(abo);

        g      = max(g, epsN);
        log_g  = log(g);

    elseif  strcmp('ded', tok{1})

        fs = 60;

        below = bif.durations_s <  out.tndt;
        bet   = bif.durations_s >= out.tndt & bif.durations_s <= points.censoring;
        abo   = bif.durations_s >  points.censoring;

        [pdf, cdf] = pdf_cdf({'ed'});

        f = @(ts, inds) out.pmix(inds)' .* pdf.ed(ts, out.theta1(inds), out.mu1(inds, :), out.tndt(inds), fs) + ...
            (1 - out.pmix(inds))' .* pdf.ed(ts, out.theta2(inds), out.mu2(inds, :), out.tndt(inds), fs);
        F = @(ts, inds) out.pmix(inds)' .* cdf.ed(ts, out.theta1(inds), out.mu1(inds, :), out.tndt(inds), fs) + ...
            (1 - out.pmix(inds))' .* cdf.ed(ts, out.theta2(inds), out.mu2(inds, :), out.tndt(inds), fs);

        epsN = 1e-12;

        if ~isempty(points.truncation)
            trunc_factor = @(inds) max(F(points.truncation, inds), epsN) ;
        else
            trunc_factor = @(inds) ones(size(ts(inds)))';
        end

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = F(points.censoring, abo) ./ trunc_factor(abo);

        g      = max(g, 1e-12);
        log_g  = log(g);

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
