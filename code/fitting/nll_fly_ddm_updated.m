function [nll, f, fd] = nll_fly_ddm_updated(params, bouts, points, model_num, iid, plot_flag, extra)
    if nargin < 7, extra = []; end
    persistent cache

    %% =========================
    % Cache model-invariant info
    % =========================
    % We must ensure the evaluator is compiled HERE, not in the sub-function
    if isempty(cache) || ~isfield(cache, model_num)
        c = struct;
        c.model_func = str2func(model_num);
        tok = regexp(model_num, '_(\D+)\d+', 'tokens');
        c.model_type = tok{1};
        
        % Instantiate model
        model = c.model_func();
        [gt, lbl] = get_ground_truth_vector(model);
        c.lbl = lbl(~isnan(gt));
        c.has_tndt = isfield(model, 'tndt');
        
        % Handles
        [pdf, cdf] = pdf_cdf({'ddm'});
        c.pdf_ddm = pdf.ddm;
        c.cdf_ddm = cdf.ddm;

        % --- COMPILE EVALUATOR ONCE ---
        % Note: This assumes 'compile_evaluate_model' does not depend on 
        % the specific length of 'bouts' but rather the field names.
        % We create a dummy structure based on input to compile the map.
        y_dummy.sm = bouts.sm; 
        y_dummy.smp = bouts.smp; 
        y_dummy.fs = bouts.fs;
        y_dummy.ln = bouts.ln; 
        y_dummy.ls = bouts.ls; 
        y_dummy.intercept = bouts.intercept;
        
        if isfield(extra,'tndt') && c.has_tndt
            model = rmfield(model,'tndt');
        end
        c.eval_fun = compile_evaluate_model(model, y_dummy, c.lbl);
        % ------------------------------

        cache.(model_num) = c;
    end
    
    c = cache.(model_num);

    %% =========================
    % Likelihood computation
    % =========================
    % If parameters are shared across flies, 'dep' and 'iid' are mathematically 
    % identical but 'dep' is slower. Only use loop if memory is an issue.
    if strcmp(iid, 'dep') && isfield(bouts, 'fly')
        flies = unique(bouts.fly);
        nll = 0;
        for i = 1:length(flies)
            mask = bouts.fly == flies(i);
            % Pass sliced data
            bouts_slice = bouts(mask,:); 
            % Note: If points are not fly-specific, this might need adjustment
            log_g = comp_loglikelihood_dddm_fast(params, bouts_slice, points, c);
            nll = nll - sum(log_g);
        end
    else
        % Default to vectorized
        log_g = comp_loglikelihood_dddm_fast(params, bouts, points, c);
        nll = -sum(log_g);
    end

    f  = [];
    fd = [];
end

function log_g = comp_loglikelihood_dddm_fast(x, bif, points, c)
    % Inputs
    ts = bif.durations_s;
    C  = points.censoring;
    t0 = points.truncation;
    epsN = 1e-12;

    % Evaluate model (Using pre-compiled c.eval_fun)
    % NOTE: Ensure c.eval_fun can handle the current size of 'bif'
    out = c.eval_fun(x); 
    
    if ~isfield(out,'tndt')
        out.tndt = zeros(size(ts));
    end

    % Masks
    bet = ts >= out.tndt & ts <= C;
    abo = ts > C;
    g = zeros(size(ts));

    % BETWEEN: pdf / truncation
    if any(bet)
        p   = out.pmix(bet);
        mu1 = out.mu1(bet); th1 = out.theta1(bet);
        mu2 = out.mu2(bet); th2 = out.theta2(bet);
        ndt = out.tndt(bet);
        t   = ts(bet);
        
        % Mixture PDF
        f1 = guard_ddm_fast(c.pdf_ddm, t, mu1, th1, ndt);
        f2 = guard_ddm_fast(c.pdf_ddm, t, mu2, th2, ndt);
        pdf_val  = p .* f1 + (1-p) .* f2;

        % Truncation Normalization
        if ~isempty(t0)
            F1 = guard_ddm_fast(c.cdf_ddm, t0, mu1, th1, ndt);
            F2 = guard_ddm_fast(c.cdf_ddm, t0, mu2, th2, ndt);
            trunc = max(1 - (p .* F1 + (1-p) .* F2), epsN);
        else
            trunc = 1;
        end
        g(bet) = pdf_val ./ trunc;
    end

    % ABOVE: survival
    if any(abo)
        p   = out.pmix(abo);
        mu1 = out.mu1(abo); th1 = out.theta1(abo);
        mu2 = out.mu2(abo); th2 = out.theta2(abo);
        ndt = out.tndt(abo);
        
        % Survival Function (1 - CDF at Censoring point)
        F1 = guard_ddm_fast(c.cdf_ddm, C, mu1, th1, ndt);
        F2 = guard_ddm_fast(c.cdf_ddm, C, mu2, th2, ndt);
        S  = 1 - (p .* F1 + (1-p) .* F2);

        if ~isempty(t0)
            Ft1 = guard_ddm_fast(c.cdf_ddm, t0, mu1, th1, ndt);
            Ft2 = guard_ddm_fast(c.cdf_ddm, t0, mu2, th2, ndt);
            trunc = max(1 - (p .* Ft1 + (1-p) .* Ft2), epsN);
        else
            trunc = 1;
        end
        g(abo) = S ./ trunc;
    end

    g = max(g, epsN);
    log_g = log(g);
end

function y = guard_ddm_fast(fun, t, mu, th, ndt)

y = zeros(size(mu));

mask = t > ndt;
if any(mask)
    if isscalar(t)
        y(mask) = fun(t, mu(mask), th(mask), ndt(mask));
    else
        y(mask) = fun(t(mask), mu(mask), th(mask), ndt(mask));
    end
end
end