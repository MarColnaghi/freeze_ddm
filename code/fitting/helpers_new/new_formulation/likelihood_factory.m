function nll_fun = likelihood_factory(model_func, surrogate, points)
% LIKELIHOOD_FACTORY  Builds a fast negative log-likelihood closure.
%
%   nll_fun = likelihood_factory(model_func, surrogate, points)
%
%   INPUTS
%     model_func  : zero-arg function that returns the model struct
%     surrogate   : table  — rows = observations, cols = predictors
%                   must have a column 'durations_s'
%     points      : struct with fields
%                     .truncation  (scalar left-truncation time t0)
%                     .censoring   (scalar right-censoring time C)
%
%   OUTPUT
%     nll_fun     : @(x) scalar NLL, x is a row-vector of coefficients

% =========================================================================
% PART A — ONE-TIME SETUP  (runs once, captured by closure)
% =========================================================================
model = model_func();
ts    = surrogate.durations_s(:);   % N x 1, force column
t0    = points.truncation;          % scalar
C     = points.censoring;           % scalar
epsN  = 1e-12;

% Flat coefficient labels (same order as x vector)
[~, lbl]  = get_ground_truth_vector_new(model);
y_vars    = surrogate.Properties.VariableNames;
y_matrix  = table2array(surrogate);     % N x P numeric matrix
N         = numel(ts);

% ------------------------------------------------------------------
% BUILD MAP : one entry per model parameter
% map.params(i).name   — field name stored in `out` struct
% map.params(i).link   — inverse-link function handle
% map.params(i).x_idx  — indices into flat x vector (coefficients)
% map.params(i).y_idx  — column indices into y_matrix  (predictors)
% ------------------------------------------------------------------
map = struct('params', []);
categories = fieldnames(model);
cnt = 1;

for i = 1:numel(categories)
    cat = categories{i};
    if ~isstruct(model.(cat)), continue; end

    params = fieldnames(model.(cat));
    for j = 1:numel(params)
        p = params{j};
        if strcmp(p, 'type') || ~isstruct(model.(cat).(p)), continue; end

        % Construct the field name that will live in `out`
        if startsWith(cat, 'process')
            eval_name = [cat '_' p];   % e.g. 'process1_mu'
        else
            eval_name = p;             % e.g. 'pmix', 'tndt'
        end

        map.params(cnt).name  = eval_name;
        map.params(cnt).link  = model.(cat).(p).link;

        block = model.(cat).(p);
        map.params(cnt).x_idx = zeros(1, numel(block.predictors));
        map.params(cnt).y_idx = zeros(1, numel(block.predictors));

        for k = 1:numel(block.predictors)
            p_name = block.predictors{k}.name;

            xi = find(strcmp(lbl,    [eval_name '_' p_name]));
            yi = find(strcmp(y_vars,  p_name));

            if isempty(xi)
                error('likelihood_factory: coefficient "%s" not found in lbl.\nAvailable labels:\n%s', ...
                    [eval_name '_' p_name], strjoin(lbl, ', '));
            end
            if isempty(yi)
                error('likelihood_factory: predictor "%s" not found in surrogate columns.\nAvailable columns:\n%s', ...
                    p_name, strjoin(y_vars, ', '));
            end

            map.params(cnt).x_idx(k) = xi;
            map.params(cnt).y_idx(k) = yi;
        end
        cnt = cnt + 1;
    end
end

% ------------------------------------------------------------------
% Sorted process names  — CRITICAL: sort so process1 < process2 < …
% and w{k} is assigned to proc_names{k} consistently.
% ------------------------------------------------------------------
proc_names = sort(categories(startsWith(categories, 'process')));
num_procs  = numel(proc_names);

% Diagnostic: print the parameter map once so you can verify field names
fprintf('\n=== likelihood_factory: parameter map ===\n');
for i = 1:numel(map.params)
    fprintf('  out.%-25s  x_idx=[%s]  y_idx=[%s]\n', ...
        map.params(i).name, ...
        num2str(map.params(i).x_idx), ...
        num2str(map.params(i).y_idx));
end
fprintf('  processes (in weight order): %s\n', strjoin(proc_names, ', '));
fprintf('=========================================\n\n');

% =========================================================================
% PART B — FAST INNER LOOP  (called thousands of times by optimizer)
% =========================================================================
nll_fun = @(x) compute_nll(x(:).');   % always pass as row vector

    function nll = compute_nll(x)

        % --- 1. Evaluate all parameters via linear predictor + link ------
        out = struct();
        for ii = 1:numel(map.params)
            m           = map.params(ii);
            % [N x K] * [K x 1] = [N x 1]
            eta         = y_matrix(:, m.y_idx) * x(m.x_idx).';
            out.(m.name) = m.link(eta(:));          % always N x 1
        end

        ndt = out.tndt(:);                          % N x 1

        % --- 2. Mixture weights ------------------------------------------
        % w{k} is N x 1 for each process k
        if num_procs == 1
            w = {ones(N, 1)};

        elseif num_procs == 2
            p_val = min(max(out.pmix(:), 1e-6), 1 - 1e-6);  % N x 1
            w     = {p_val, 1 - p_val};

        else
            % For 3+ processes use softmax on pmix_1 … pmix_(K-1)
            % (Extend here if your model ever has >2 processes.)
            error('likelihood_factory: num_procs > 2 not yet implemented.');
        end

        % --- 3. Aggregate pdf/cdf across processes -----------------------
        % ts is N x 1 ; t0, C are scalars
        % dist_dispatch must return N x 1 for any time argument because
        % ndt is N x 1 (observation-specific).

        f_mix  = zeros(N, 1);
        F_mix  = zeros(N, 1);   %#ok<NASGU>  (kept for possible future use)
        F0_mix = zeros(N, 1);   % CDF at truncation point t0
        FC_mix = zeros(N, 1);   % CDF at censoring point C

        for k = 1:num_procs
            p_name = proc_names{k};
            type   = model.(p_name).type;

            % Observed times
            [fk, Fk]   = dist_dispatch(type, out, p_name, ndt, ts);

            % Truncation point (scalar t0  →  N x 1 because ndt is N x 1)
            [~,  F0k]  = dist_dispatch(type, out, p_name, ndt, t0);

            % Censoring point (scalar C)
            [~,  FCk]  = dist_dispatch(type, out, p_name, ndt, C);

            wk = w{k};                              % N x 1

            f_mix  = f_mix  + wk .* fk(:);
            F_mix  = F_mix  + wk .* Fk(:);         %#ok<AGROW>
            F0_mix = F0_mix + wk .* F0k(:);
            FC_mix = FC_mix + wk .* FCk(:);
        end

        % --- 4. Truncation factor ----------------------------------------
        trunc_factor = max(1 - F0_mix, epsN);      % N x 1

        % --- 5. Likelihood per observation --------------------------------
        g     = zeros(N, 1);
        valid = (ts >= ndt + 0.001);                % below-ndt guard
        bet   = valid & (ts <= C);                  % observed (between)
        abo   = (ts >  C);                          % right-censored (above)

        g(bet) = f_mix(bet)        ./ trunc_factor(bet);
        g(abo) = (1 - FC_mix(abo)) ./ trunc_factor(abo);

        % --- 6. Guard underflow & sum ------------------------------------
        g   = max(g, epsN);
        nll = -sum(log(g));

        % --- 7. Soft penalty: ndt must stay below min observed RT --------
        violators = ndt >= min(ts);
        if any(violators)
            nll = nll + 1e6 * sum(violators);
        end

    end % compute_nll

end % likelihood_factory