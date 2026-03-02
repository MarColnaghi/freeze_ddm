function [lb, plb, pub, ub, prior_info] = extract_bounds_from_model_new(model)
    
% Initialize as empty row vectors
    lb = [];
    ub = [];
    % Initialize prior_info as a struct array with the correct fields
    prior_info = struct('type', {}, 'lambda', {});
    
    categories = fieldnames(model);
    
    for i = 1:numel(categories)
        catName = categories{i};
        category = model.(catName);
        
        % Check if this category is actually a struct (e.g., process1, mixing)
        if ~isstruct(category), continue; end
        
        params = fieldnames(category);
        for j = 1:numel(params)
            pName = params{j};
            
            % Skip 'type' string or other non-struct metadata
            if strcmp(pName, 'type') || ~isstruct(category.(pName))
                continue; 
            end
            
            block = category.(pName);
            
            % Ensure 'predictors' exists and is a cell array
            if ~isfield(block, 'predictors') || ~iscell(block.predictors)
                continue;
            end
            
            numPreds = numel(block.predictors);
            block_lb = zeros(1, numPreds);
            block_ub = zeros(1, numPreds);
            
            % Build block_prior individually to avoid concatenation mismatches
            for k = 1:numPreds
                pred = block.predictors{k};
                
                % 1. Extract Bounds
                block_lb(k) = pred.bounds(1);
                block_ub(k) = pred.bounds(2);
                
                % 2. Extract Prior (with defaults)
                this_prior = struct('type', 'trapezoidal', 'lambda', []);
                if isfield(pred, 'prior')
                    this_prior.type = pred.prior;
                    if isfield(pred, 'lambda')
                        this_prior.lambda = pred.lambda;
                    end
                end
                
                % Append prior_info one by one for safety
                prior_info(end+1) = this_prior; 
            end
            
            % Append bounds
            lb = [lb, block_lb];
            ub = [ub, block_ub];
        end
    end

    % --- Generate Plausible Bounds ---
    range = ub - lb;
    % Shift inwards by 10%
    shift = range * 0.1;
    
    % Ensure a minimum shift for non-fixed parameters
    shift(range > 0) = max(shift(range > 0), 0.01);
    
    plb = lb + shift;
    pub = ub - shift;
    
    % Handle fixed parameters (LB == UB)
    fixed_idx = (lb == ub);
    plb(fixed_idx) = lb(fixed_idx);
    pub(fixed_idx) = ub(fixed_idx);
end