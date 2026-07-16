function [pdf, surv, t, x, P] = fpe_cn(drift, prm)
% FPE_CN  Crank-Nicolson Fokker-Planck solver for a leaky accumulator (teaching
%         implementation, pure MATLAB -- the readable twin of
%         fitting/instantaneous/leaky_pde_robust.cpp).
%
%   [pdf, surv, t, x] = fpe_cn(drift, prm)
%
%   drift : per-frame drift RATE mu(t) (vector, length T).
%   prm   : struct with dt, dx, sigma_sq, x0, x_min, bound, lambda.
%
%   pdf  : first-passage density on t = (0:T-1)*dt   (pdf(1) = 0, see step 5)
%   surv : probability still un-crossed at t = T*dt  (censoring mass)
%   t, x : the time and space grids.
%   P    : OPTIONAL [T x n] density p(x,t) at every step -- only stored when
%          asked for (it is T*n doubles). Use it to SEE the method work:
%          sims/demo_fpe_anatomy.m plots it.
%
% ============================ THE DERIVATION ============================
%
% 1. SDE -> FOKKER-PLANCK.  The accumulator
%
%        dx = (mu(t) - lambda*x) dt + sigma dW
%
%    has a density p(x,t) obeying the Fokker-Planck (forward Kolmogorov) equation
%
%        dp/dt = -d/dx[ v(x,t) p ]  +  (sigma^2/2) d2p/dx2,    v = mu(t) - lambda*x
%                \_____________/       \________________/
%                  advection             diffusion
%
%    Note v is STATE-DEPENDENT: the leak term -lambda*x is why v is a vector over
%    the grid, not a scalar. lambda = 0 makes v spatially constant.
%
% 2. DISCRETISE IN SPACE (centred differences on a grid x_i, spacing dx):
%
%        d2p/dx2 |_i  ~  ( p_{i-1} - 2 p_i + p_{i+1} ) / dx^2
%        d/dx(v p)|_i ~  ( v_{i+1} p_{i+1} - v_{i-1} p_{i-1} ) / (2 dx)
%
% 3. CRANK-NICOLSON IN TIME.  Average the right-hand side over steps n and n+1
%    (half explicit, half implicit). That "half" is where every /2 below comes
%    from, and it buys 2nd-order accuracy in dt AND unconditional stability:
%
%        (p^{n+1}_i - p^n_i)/dt = 1/2 [ RHS^{n+1} + RHS^n ]_i
%
%    Multiplying out and collecting terms gives exactly two constants:
%
%        C    = (sigma^2/2) * (dt/2) / dx^2   =  sigma_sq*dt / (4*dx^2)
%        beta_i = v_i * (dt/2) / (2*dx)       =  v_i * dt / (4*dx)
%
%    ...which is where the mysterious /4 in both comes from. (Same constants
%    appear as C_diff and drift_mult in leaky_pde_robust.cpp:62-63, and as
%    alpha / drift_mult in bayes_fpe tv_core.py:107,101.)
%
% 4. THE TRIDIAGONAL SYSTEM.  Implicit terms to the left, explicit to the right:
%
%      A p^{n+1} = B p^n, with row i of each:
%        A:  (-C - beta_{i-1}) p_{i-1} + (1 + 2C) p_i + (-C + beta_{i+1}) p_{i+1}
%        B:  ( C + beta_{i-1}) p_{i-1} + (1 - 2C) p_i + ( C - beta_{i+1}) p_{i+1}
%
%    A is tridiagonal, so each step costs O(n) via the Thomas algorithm.
%
%    BOUNDARIES:
%      upper (x = bound): ABSORBING   -> p_n = 0. This is the decision threshold.
%      lower (x = x_min): reflecting, imposed here as zero-gradient p_1 = p_2.
%        (Only a containment wall, NOT a second decision bound. NB bayes_fpe uses
%        a proper flux-balance row instead; zero-gradient zeroes the diffusive
%        flux only, so the two differ slightly when v ~= 0 at x_min.)
%
% 5. THE FIRST-PASSAGE DENSITY.  The survival function is S(t) = int p(x,t) dx.
%    The lower edge reflects, so the ONLY way mass leaves is by being absorbed at
%    the threshold. Therefore
%
%        f(t) = -dS/dt  ~  ( S(t) - S(t+dt) ) / dt
%
%    i.e. the mass that vanished this step, per unit time. That is the whole
%    trick, and it is what line "pdf(k) = ..." below computes. (The textbook
%    alternative is the probability flux J = v p - (sigma^2/2) dp/dx evaluated at
%    x = bound; mass-loss is equivalent and better behaved numerically.)
%    pdf(1) = 0 because no mass can be absorbed before the first step.
%
% CONVENTION: step k advances t=(k-2)*dt -> t=(k-1)*dt and uses drift(k), i.e.
% the drift at the END of the step. This matches leaky_pde_robust.cpp exactly
% (it reads drifts[k] 0-based for the step ending at k*dt). sim_leaky_accumulator
% instead uses the drift at the START of the step -- a one-frame offset on a
% time-varying drive. See sims/test_sim_vs_fpe.m.
% =======================================================================

    dt = prm.dt;  dx = prm.dx;  sigma_sq = prm.sigma_sq;  lambda = prm.lambda;

    % ---- grid --------------------------------------------------------------
    n = round((prm.bound - prm.x_min)/dx) + 1;
    x = (prm.x_min + (0:n-1)'*dx);
    T = numel(drift);
    t = (0:T-1)' * dt;

    % ---- the two CN constants (step 3) -------------------------------------
    C  = sigma_sq * dt / (4 * dx^2);
    dm = dt / (4 * dx);
    v_leak = -lambda * x;              % state-dependent part of v; scalar mu adds on

    % ---- initial condition -------------------------------------------------
    % Default: all mass at x0 as a DELTA (height 1/dx so it integrates to 1),
    % snapped to the nearest node. This matches leaky_pde_robust.cpp exactly.
    %
    % CAVEAT worth knowing: CN does not damp high-frequency modes (its
    % amplification factor tends to -1, not 0), so a delta start RINGS -- the
    % density develops an alternating-sign checkerboard around x0 that never
    % dies. It is harmless for the FPT density (the wiggles are symmetric, cancel
    % in the mass integral, and sit far from the bound -- KS vs the exact Wald is
    % still ~1e-4), but it is ugly, and it is why you must NOT clip p to >= 0:
    % clipping the negative lobes creates mass.
    %
    % Set prm.init_sigma_dx = w to start from a GAUSSIAN bump of width w*dx
    % instead. Smooth initial data has little high-frequency content, so it does
    % not ring. This is exactly what bayes_fpe does (tv_core.py:77-86), and why
    % bayes_fpe CAN afford to clip.
    start_idx = round((prm.x0 - prm.x_min)/dx) + 1;
    start_idx = min(max(start_idx, 1), n);
    if isfield(prm, 'init_sigma_dx') && ~isempty(prm.init_sigma_dx)
        s0 = max(prm.init_sigma_dx * dx, 1e-12);
        p  = exp(-0.5 * ((x - x(start_idx)) / s0).^2);
        p(end) = 0;                    % the absorbing node holds no mass
        p  = p / max(sum(p)*dx, 1e-30);
    else
        p = zeros(n,1);
        p(start_idx) = 1/dx;
    end

    pdf   = zeros(T,1);                % pdf(1) = 0 (nothing absorbed yet)
    sum_p = sum(p);

    keep_P = nargout > 4;              % only pay the memory when asked
    if keep_P, P = zeros(T, n); P(1,:) = p.'; else, P = []; end

    n_rann = 0;
    if isfield(prm, 'rannacher') && ~isempty(prm.rannacher), n_rann = prm.rannacher; end

    for k = 2:T
        v    = drift(k) + v_leak;      % v_i = mu(t) - lambda*x_i   (end-of-step mu)

        % ---- RANNACHER STARTUP -------------------------------------------
        % CN is A-stable but NOT L-stable: its amplification factor for the
        % highest modes tends to -1, so they flip sign every step and never
        % decay. A non-smooth start (a delta) is pure high-frequency content, so
        % it RINGS -- the checkerboard around x0. Refining dx makes it worse, not
        % better, because C = sigma^2*dt/(4dx^2) grows.
        % Backward Euler IS L-stable (g = 1/(1+2C(1-cos)) -> 0), so it
        % annihilates those modes. Taking the first few steps with BE and then
        % switching to CN kills the ringing at its source while keeping CN's
        % O(dt^2) accuracy overall (a couple of first-order steps contribute only
        % O(dt^2) globally). Rannacher (1984); Giles & Carter (2006).
        % BE uses the FULL dt where CN uses dt/2, hence the doubling.
        if k-1 <= n_rann
            Ck = 2*C;  betak = v * (2*dm);  is_be = true;
        else
            Ck = C;    betak = v * dm;      is_be = false;
        end
        beta = betak;

        % ---- assemble A (implicit) and the RHS, step 4 ---------------------
        % interior rows 2..n-1
        lo = -Ck - beta(1:n-2);        % sub-diagonal: multiplies p_{i-1}
        mid = (1 + 2*Ck) * ones(n,1);  % main diagonal
        up = -Ck + beta(3:n);          % super-diagonal: multiplies p_{i+1}

        rhs = zeros(n,1);
        if is_be
            % backward Euler: (I - dt*L) p^{n+1} = p^n  -- RHS is just p^n
            rhs(2:n-1) = p(2:n-1);
        else
            % Crank-Nicolson: RHS is B*p^n
            rhs(2:n-1) = ( C + beta(1:n-2)) .* p(1:n-2) ...
                       + (1 - 2*C)          .* p(2:n-1) ...
                       + ( C - beta(3:n))   .* p(3:n);
        end

        % ---- boundary rows -------------------------------------------------
        % lower: zero-gradient  ->  p_1 - p_2 = 0
        mid(1) = 1;  up_1 = -1;  rhs(1) = 0;
        % upper: absorbing      ->  p_n = 0
        mid(n) = 1;  lo_n = 0;   rhs(n) = 0;

        % pack the diagonals for the Thomas solve: row i is
        %   a(i-1)*p_{i-1} + b(i)*p_i + c(i)*p_{i+1} = d(i)
        a = [lo; lo_n];                % length n-1, a(i) belongs to row i+1
        b = mid;                       % length n
        c = [up_1; up];                % length n-1, c(i) belongs to row i
        p_new = thomas(a, b, c, rhs);
        % NB do NOT clip p_new to >= 0 here. CN ringing around the delta initial
        % condition makes some nodes slightly negative; those negatives are
        % symmetric and cancel, so the total mass still telescopes exactly.
        % Clipping them to zero CREATES mass and destroys the identity
        % int(pdf)dt + surv = 1 (it silently inflated the total to ~1.97).
        % leaky_pde_robust.cpp does not clip either. (bayes_fpe can afford to,
        % because it starts from a smooth Gaussian bump, not a delta.)

        % ---- step 5: absorbed mass this step IS the density ----------------
        sum_p_new = sum(p_new);
        pdf(k)    = max(0, (sum_p - sum_p_new)) * dx / dt;

        sum_p = sum_p_new;
        p     = p_new;
        if keep_P, P(k,:) = p.'; end
    end

    surv = sum(p) * dx;                % mass never absorbed within the window
end

% =========================================================================
function x = thomas(a, b, c, d)
% Thomas algorithm (tridiagonal Gaussian elimination), O(n).
%   row i:  a(i-1)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i)
% Mirrors solve_tdma() in leaky_pde_robust.cpp.
    n  = numel(b);
    cp = zeros(n-1,1);
    dp = zeros(n,1);

    cp(1) = c(1) / b(1);               % forward sweep: eliminate the sub-diagonal
    dp(1) = d(1) / b(1);
    for i = 2:n-1
        m     = 1 / (b(i) - a(i-1)*cp(i-1));
        cp(i) = c(i) * m;
        dp(i) = (d(i) - a(i-1)*dp(i-1)) * m;
    end
    dp(n) = (d(n) - a(n-1)*dp(n-1)) / (b(n) - a(n-1)*cp(n-1));

    x    = zeros(n,1);                 % back substitution
    x(n) = dp(n);
    for i = n-1:-1:1
        x(i) = dp(i) - cp(i)*x(i+1);
    end
end
