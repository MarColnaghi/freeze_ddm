function out = evaluate_model_new(model, x, y)
% ------------------------------------------------------------
% Custom Naming Logic:
% - process1/2 -> "process1_mu_sm"
% - mixing     -> "pmix_smp"
% - shared     -> "tndt_intercept"
% ------------------------------------------------------------
    y_vars = y.Properties.VariableNames;
    x_vars = x.Properties.VariableNames;
    N = height(y);
    out = table();
    
    categories = fieldnames(model);
    
    for i = 1:numel(categories)
        catName = categories{i};
        category = model.(catName);
        params = fieldnames(category);
        
        for j = 1:numel(params)
            pName = params{j};
            
            % Skip metadata fields
            if strcmp(pName, 'type') || ~isstruct(category.(pName))
                continue; 
            end
            
            block = category.(pName);
            
            % --- CUSTOM NAMING LOGIC ---
            % If it's a process, keep the category prefix. 
            % If it's mixing/shared, use only the parameter name.
            if startsWith(catName, 'process')
                fullParamName = [catName '_' pName];
            else
                fullParamName = pName;
            end
            % ---------------------------

            eta = zeros(N, 1);
            
            for k = 1:numel(block.predictors)
                pred_name = block.predictors{k}.name;
                
                % 1. Get data from y
                Xcol = y.(pred_name);
                
                % 2. Look up coefficient in x (e.g., "process1_mu_sm" or "pmix_smp")
                coeff_name = [fullParamName '_' pred_name];
                
                if ~ismember(coeff_name, x_vars)
                    error('Coefficient "%s" not found in x.', coeff_name);
                end
                beta = x.(coeff_name); 

                % ---- Handle scalar vs time-series predictors ----
                if iscell(Xcol)
                    Xi = cell2mat(Xcol); 
                    if size(beta, 2) == 1 
                        eta = eta + Xi * beta;
                    else 
                        eta = eta + sum(Xi .* beta, 2); 
                    end
                else
                    eta = eta + Xcol .* beta;
                end
            end
            
            % Apply link and store in output table
            out.(fullParamName) = block.link(eta);
        end
    end
end