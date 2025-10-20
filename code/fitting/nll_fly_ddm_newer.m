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
    tbl.durations_s = [0:0.005:62]';
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

        g = comp_loglikelihood(params, tbl, points, model_func, iid, tok, extra);
        %[~, G] = comp_loglikelihood(params, tbl, points, model_func, iid, tok, extra);

        F = F + G;
        f = f + exp(g);

    end

    f = f ./ height(bouts);
    F = F ./ height(bouts);

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


    elseif  strcmp('dddm', tok{1})

        bet = bif.durations_s > out.tndt & bif.durations_s <= points.censoring;

        % pdf and cdf for single bound ddm
        [pdf, cdf] = pdf_cdf({'ddm', 'exp'});

        f = @(ts, inds) out.pmix(inds) .* pdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds)) + ...
            (1 - out.pmix(inds)) .* pdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));
        F = @(ts, inds) out.pmix(inds) .* cdf.ddm(ts, out.mu1(inds), out.theta1(inds), out.tndt(inds)) + ...
            (1 - out.pmix(inds)) .* cdf.ddm(ts, out.mu2(inds), out.theta2(inds), out.tndt(inds));

        if ~isempty(points.truncation) && points.truncation > out.tndt(1)
            trunc_factor = @(inds) 1 - F(points.truncation, inds);
        else
            trunc_factor = @(inds) ones(size(ts(inds)));
        end

        g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
        g(abo) = (1 - F(points.censoring, abo)) ./ trunc_factor(abo);
        

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
    end

end

g = max(g, 1e-6);
log_g = log(g); %- log(trunc_factor(ones(size(ts))));

end

