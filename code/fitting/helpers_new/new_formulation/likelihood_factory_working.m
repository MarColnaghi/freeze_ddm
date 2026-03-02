function nll_fun = likelihood_factory_working(model_func, surrogate, points)
% --- PART A: HIERARCHICAL SETUP (Runs Once) ---
model = model_func();
ts = surrogate.durations_s(:);
t0 = points.truncation;
C  = points.censoring;
epsN = 1e-12;

y_vars = surrogate.Properties.VariableNames;
y_matrix = table2array(surrogate);
[~, lbl] = get_ground_truth_vector_new(model);

% Identify processes (e.g., process1, process2)
categories = fieldnames(model);
proc_names = categories(startsWith(categories, 'process'));
num_procs = numel(proc_names);

% --- 1. Map Mixture Weights (pmix) ---
pmix_info = [];
if num_procs > 1
    if isfield(model, 'shared') && isfield(model.shared, 'pmix')
        % Shared parameters usually don't have the "shared_" prefix in lbl
        pmix_info = map_parameter('', 'pmix', model.shared.pmix, lbl, y_vars);
    elseif isfield(model, 'pmix') && isfield(model.pmix, 'pmix')
        % Alternative structure check
        pmix_info = map_parameter('pmix', 'pmix', model.pmix.pmix, lbl, y_vars);
    end
end

% --- 2. Build Process Recipes ---
recipes = struct();
for k = 1:num_procs
    p_name = proc_names{k};
    dist_type = model.(p_name).type;
    [pdf_set, cdf_set] = pdf_cdf({dist_type});
    
    recipes(k).pdf_fh = pdf_set.(dist_type);
    recipes(k).cdf_fh = cdf_set.(dist_type);
    
    % Get functional parameters (exclude metadata and timing)
    fields = fieldnames(model.(p_name));
    dist_params = fields(~ismember(fields, {'type', 'tndt'}));
    
    % Map standard parameters (mu, theta, lambda, etc.)
    for p = 1:numel(dist_params)
        p_key = dist_params{p};
        recipes(k).params(p) = map_parameter(p_name, p_key, model.(p_name).(p_key), lbl, y_vars);
    end
    
    % --- 3. Handle NDT (Local in Process or Shared) ---
    if isfield(model.(p_name), 'tndt')
        % NDT is inside the process (model.process1.tndt)
        recipes(k).ndt_info = map_parameter(p_name, 'tndt', model.(p_name).tndt, lbl, y_vars);
    elseif isfield(model, 'shared') && isfield(model.shared, 'tndt')
        % NDT is shared (model.shared.tndt)
        recipes(k).ndt_info = map_parameter('', 'tndt', model.shared.tndt, lbl, y_vars);
    elseif isfield(model, 'tndt') && isfield(model.tndt, 'tndt')
        % NDT is in its own top-level block
        recipes(k).ndt_info = map_parameter('tndt', 'tndt', model.tndt.tndt, lbl, y_vars);
    else
        error('Non-decision time (tndt) mapping not found for %s.', p_name);
    end
end

% --- PART B: THE AGNOSTIC INNER LOOP ---
nll_fun = @(x) compute_nll(x(:)');

    function nll = compute_nll(x)
        % 1. Evaluate Mixture Weights
        if num_procs > 1 && ~isempty(pmix_info)
            pm = pmix_info.link(y_matrix(:, pmix_info.y_idx) * x(pmix_info.x_idx).');
            pm = min(max(pm(:), 1e-6), 1 - 1e-6);
            w = {pm, 1 - pm};
        else
            % Default weight for single process or unweighted models
            w = repmat({ones(size(ts))}, 1, num_procs);
        end

        f_mix = zeros(size(ts)); F_mix = zeros(size(ts));
        F0_mix = zeros(size(ts)); FC_mix = zeros(size(ts));
        
        % Track NDT to determine likelihood validity
        min_ndt = inf(size(ts)); 

        for k = 1:num_procs
            % Extract Distribution Parameters
            num_p = numel(recipes(k).params);
            p_vals = cell(1, num_p);
            for p = 1:num_p
                eta = y_matrix(:, recipes(k).params(p).y_idx) * x(recipes(k).params(p).x_idx).';
                p_vals{p} = recipes(k).params(p).link(eta);
            end
            
            % Extract NDT for this component
            ndt_k = recipes(k).ndt_info.link(y_matrix(:, recipes(k).ndt_info.y_idx) * x(recipes(k).ndt_info.x_idx).');
            min_ndt = min(min_ndt, ndt_k); 
            
            % Compute Likelihoods
            fk  = guard_generic(recipes(k).pdf_fh, ts, p_vals{:}, ndt_k);
            Fk  = guard_generic(recipes(k).cdf_fh, ts, p_vals{:}, ndt_k);
            F0k = guard_generic(recipes(k).cdf_fh, t0, p_vals{:}, ndt_k);
            FCk = guard_generic(recipes(k).cdf_fh, C,  p_vals{:}, ndt_k);
            
            % Sum to mixture
            f_mix  = f_mix  + w{k} .* fk(:);
            F_mix  = F_mix  + w{k} .* Fk(:);
            F0_mix = F0_mix + w{k} .* F0k(:);
            FC_mix = FC_mix + w{k} .* FCk(:);
        end
        
        % 3. Global Truncation & Censoring
        trunc = max(1 - F0_mix, epsN);
        g = zeros(size(ts));
        
        % Logical indexing for speed and boundary safety
        valid = (ts > min_ndt + 0.0001); 
        bet = valid & (ts <= C);
        abo = (ts > C);
        
        g(bet) = f_mix(bet) ./ trunc(bet);
        g(abo) = (1 - FC_mix(abo)) ./ trunc(abo);
        
        % 4. Final Negative Log Likelihood
        nll = -sum(log(max(g, epsN)));
    end

    % --- Helper: Map parameters to X (optimization) and Y (data) indices ---
    function info = map_parameter(cat_name, p_name, block, lbl, y_vars)
        info.link = block.link;
        num_preds = numel(block.predictors);
        info.x_idx = zeros(1, num_preds);
        info.y_idx = zeros(1, num_preds);
        
        % Construct naming: 'process1_mu' or just 'pmix' if shared
        if isempty(cat_name)
            prefix = p_name;
        else
            prefix = [cat_name '_' p_name];
        end
        
        for m = 1:num_preds
            pred_name = block.predictors{m}.name;
            target_name = [prefix '_' pred_name];
            
            % Find in optimization vector labels
            x_match = find(strcmp(lbl, target_name));
            if isempty(x_match)
                error('Likelihood Factory: Parameter "%s" not found in model labels.', target_name);
            end
            info.x_idx(m) = x_match;
            
            % Find in data table variables
            y_match = find(strcmp(y_vars, pred_name));
            if isempty(y_match)
                error('Likelihood Factory: Predictor "%s" not found in data table.', pred_name);
            end
            info.y_idx(m) = y_match;
        end
    end

    % --- Helper: Safeguard PDF/CDF calls ---
    function val = guard_generic(fh, t, varargin)
        try
            val = fh(t, varargin{:});
            val(isnan(val) | isinf(val)) = 0;
        catch
            val = zeros(size(t));
        end
    end
end