function logL = loglik_ddm_binned(ts, bet, abo, F, points)

% Discrete-time likelihood evaluation for a continuous DDM model
% ts     : reaction times (frame-based)
% bet    : logical index of observed (not censored) trials
% abo    : logical index of right-censored trials
% F      : handle @(ts, inds) continuous CDF (mixture already inside)
% points : struct with fields
%          .dt          sampling interval
%          .truncation  left truncation time
%          .censoring   censoring time

    dt   = points.dt;
    t0   = points.truncation;
    C    = points.censoring;
    epsN = 1e-12;

    logL = 0;

    % -------- Truncation factor (discrete!) --------
    Z = 1 - F(t0 - dt, true(size(ts)));
    Z = max(Z, epsN);

    % -------- Observed trials --------
    if any(bet)
        t_prev = ts(bet) - dt;

        Fk   = F(ts(bet),     bet);
        Fk_1 = F(t_prev,      bet);

        pk = Fk - Fk_1;
        pk = max(pk, epsN);

        logL = logL + sum(log(pk)) - sum(log(Z(bet)));
    end

    % -------- Right-censored trials --------
    if any(abo)
        Sk = 1 - F(C, abo);
        Sk = max(Sk, epsN);

        logL = logL + sum(log(Sk)) - sum(log(Z(abo)));
    end
end
