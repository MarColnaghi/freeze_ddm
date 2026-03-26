function [f1, f2, fd ] = extract_individual_processes(params, bouts, points, model_num, iid, plot_flag, extra)

if nargin < 7
    extra = [];
end

model_func = str2func(model_num);
pattern = '_(\D+)\d+';
tok = regexp(model_num, pattern, 'tokens');

if strcmp(plot_flag, 'p')


    fd = (0:1/60:(points.censoring + 1/60))';
    num_times = height(fd);
    num_bouts = height(bouts);

    tbl = table();
    tbl.durations_s = fd;
    tbl.sm = zeros(num_times, 1);
    tbl.smp = zeros(num_times, 1);
    tbl.fs = zeros(num_times, 1);
    tbl.ln = zeros(num_times, 1);
    tbl.ls = zeros(num_times, 1);
    tbl.tsl = zeros(num_times, 1);
    tbl.intercept = zeros(num_times, 1);

    f1 = zeros(num_times, 1);
    f2 = zeros(num_times, 1);

    hold on

    ec = extra;

    for idx_bout = 1:num_bouts

        tbl.sm(:) = bouts.sm(idx_bout);
        tbl.smp(:) = bouts.smp(idx_bout);
        tbl.fs(:) = bouts.fs(idx_bout);
        tbl.ln(:) = bouts.ln(idx_bout);
        tbl.ls(:) = bouts.ls(idx_bout);
        tbl.tsl(:) = bouts.tsl(idx_bout);
        tbl.intercept(:) = bouts.intercept(idx_bout);

        if strcmp('ed', tok{1}) || strcmp('ded', tok{1})
            ec.soc_mot_array = extra.soc_mot_array(idx_bout, :);
        end
        [log_g1, log_g2] = comp_loglikelihood(params, tbl, points, model_func, iid, tok, ec);

        f1 = f1 + exp(log_g1);
        f2 = f2 + exp(log_g2);

    end

    f1 = f1 ./ num_bouts;
    f2 = f2 ./ num_bouts;

    f2(fd < points.truncation) = 0;
    f1(fd < points.truncation) = 0;

    %     if strcmp('ed', tok{1}) || strcmp('ded', tok{1})
    %         f(1:end-1) = f(1:end-1)* 60;
    %     end

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

function [log_g1, log_g2] = comp_loglikelihood(x, bouts_individual_fly, points, model_func, iid, tok, extra)

t0   = points.truncation;
C    = points.censoring;
epsN = 1e-12;
fs = 60;

bif = bouts_individual_fly;

ts = bif.durations_s;
y = table;
y.sm = bif.sm;
y.smp = bif.smp;
y.fs = bif.fs;
y.ln = bif.ln;
y.ls = bif.ls;
y.tsl = bif.tsl;

y.intercept = bif.intercept;

model = model_func();
if isfield(extra, 'fixed_param')
    model = rmfield(model, fieldnames(extra.fixed_param));
end

[gt, lbl] = get_ground_truth_vector(model);
lbl = lbl(~isnan(gt));
gt_table = array2table(x, 'VariableNames', lbl);

if strcmp('ed', tok{1}) || strcmp('ded', tok{1}) || strcmp('sddmtv', tok{1})
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

below = bif.durations_s < out.tndt;
bet = bif.durations_s >= out.tndt & bif.durations_s <= C;
abo = bif.durations_s > C;

g1 = zeros(size(ts));
g2 = zeros(size(ts));

if strcmp('iid', iid)

    if  strcmp('dddm', tok{1})

        [pdf, cdf] = pdf_cdf({'ddm'});

        pdf_ddm_raw = pdf.ddm;
        cdf_ddm_raw = cdf.ddm;

        pdf.ddm = @(ts, mu, theta, ndt) guard_ddm(pdf_ddm_raw, ts, mu, theta, ndt);
        cdf.ddm = @(ts, mu, theta, ndt) guard_ddm(cdf_ddm_raw, ts, mu, theta, ndt);

        f1 = @(ts, inds) out.pmix(inds) .* pdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds));
        f2 = @(ts, inds) (1 - out.pmix(inds)) .* pdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));

        F1 = @(ts, inds) out.pmix(inds) .* cdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds));
        F2 = @(ts, inds) (1 - out.pmix(inds)) .* cdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));

        F = @(ts, inds) F1(ts, inds) + F2(ts, inds);
        trunc_factor = @(inds) max(1 - F(t0, inds), epsN);

        g1(bet) = f1(ts(bet), bet) ./ trunc_factor(bet);
        g1(abo) = (out.pmix(abo) - F1(C, abo)) ./ trunc_factor(abo);

        g1 = max(g1, epsN);
        log_g1= log(g1);

        g2(bet) = f2(ts(bet), bet) ./ trunc_factor(bet);
        g2(abo) = ((1 - out.pmix(abo)) - F2(C, abo)) ./ trunc_factor(abo);

        g2 = max(g2, epsN);
        log_g2= log(g2);


    elseif  strcmp('ded', tok{1})

       fs = 60;
        % Use the same index definitions as your joint code
        below = bif.durations_s <  out.tndt;
        bet   = bif.durations_s >= out.tndt & bif.durations_s <= points.censoring;
        abo   = bif.durations_s >  points.censoring;
        
        [pdf, cdf] = pdf_cdf({'ed'});
        
        % Component 1 functions
        f1 = @(ts, inds) out.pmix(inds)' .* pdf.ed(ts, out.theta1(inds), out.mu1(inds, :), out.tndt(inds), fs);
        F1 = @(ts, inds) out.pmix(inds)' .* cdf.ed(ts, out.theta1(inds), out.mu1(inds, :), out.tndt(inds), fs);
        
        % Component 2 functions
        f2 = @(ts, inds) (1 - out.pmix(inds))' .* pdf.ed(ts, out.theta2(inds), out.mu2(inds, :), out.tndt(inds), fs);
        F2 = @(ts, inds) (1 - out.pmix(inds))' .* cdf.ed(ts, out.theta2(inds), out.mu2(inds, :), out.tndt(inds), fs);
        
        % Full Mixture (Only for the denominator)
        F_total = @(ts, inds) F1(ts, inds) + F2(ts, inds);
        
        % MATCH THE JOINT TRUNCATION LOGIC
        if ~isempty(points.truncation)
            % Use F_total here so both g1 and g2 are scaled by the same mixture mass
            trunc_denom = @(inds) max(F_total(points.truncation - 1/fs, inds), 1e-12);
        else
            trunc_denom = @(inds) ones(size(inds))';
        end
        
        % --- Component 1 ---
        g1(bet) = f1(ts(bet), bet) ./ trunc_denom(bet) * fs; % Scale by fs
        g1(abo) = F1(points.censoring, abo) ./ trunc_denom(abo); % Match joint 'F' logic
        g1 = max(g1, 1e-12);
        log_g1 = log(g1);
        
        % --- Component 2 ---
        g2(bet) = f2(ts(bet), bet) ./ trunc_denom(bet) * fs; % Scale by fs
        g2(abo) = F2(points.censoring, abo) ./ trunc_denom(abo); % Match joint 'F' logic
        g2 = max(g2, 1e-12);
        log_g2 = log(g2);

    elseif  strcmp('m', tok{1})

        [pdf, cdf] = pdf_cdf({'ddm'});

        pdf_ddm_raw = pdf.ddm;
        cdf_ddm_raw = cdf.ddm;

        pdf.ddm = @(ts, mu, theta, ndt) guard_ddm(pdf_ddm_raw, ts, mu, theta, ndt);
        cdf.ddm = @(ts, mu, theta, ndt) guard_ddm(cdf_ddm_raw, ts, mu, theta, ndt);

        f1 = @(ts, inds) out.pmix(inds) .* pdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds));
        f2 = @(ts, inds) (1 - out.pmix(inds)) .* pdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));

        F = @(ts, inds) out.pmix(inds) .* cdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds)) + ...
            (1 - out.pmix(inds)) .* cdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));

        trunc_factor = @(inds) max(1 - F(t0, inds), epsN);

        g1(bet) = f1(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = (1 - F(C, abo)) ./ trunc_factor(abo);

        g1 = max(g1, epsN);
        log_g1= log(g1);

        g2(bet) = f2(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = (1 - F(C, abo)) ./ trunc_factor(abo);

        g2 = max(g2, epsN);
        log_g2= log(g2);

    
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
