function [xtrue, lbl, mask] = get_ground_truth_vector_new(model)
% ------------------------------------------------------------
% Generates ground truth vector and labels for nested model
% Logic:
% - process1/2 -> "process1_mu_sm"
% - mixing/shared -> "pmix_smp"
% ------------------------------------------------------------
    categories = fieldnames(model);
    
    xtrue_parts = {};
    lbl_parts = {};
    mask_parts = {};

    for i = 1:numel(categories)
        catName = categories{i};
        category = model.(catName);
        params = fieldnames(category);
        
        for j = 1:numel(params)
            pName = params{j};
            
            % Skip metadata (like 'type')
            if strcmp(pName, 'type') || ~isstruct(category.(pName))
                continue; 
            end
            
            block = category.(pName);
            
            % --- CUSTOM NAMING LOGIC (Match evaluate_model) ---
            if startsWith(catName, 'process')
                namePrefix = [catName '_' pName];
            else
                namePrefix = pName;
            end
            
            % Iterate through predictors in this block
            numPreds = numel(block.predictors);
            block_xtrue = zeros(1, numPreds);
            block_lbl = cell(1, numPreds);
            
            for k = 1:numPreds
                predName = block.predictors{k}.name;
                
                % Store ground truth value
                block_xtrue(k) = block.ground_truth(k);
                
                % Create label: e.g., "process1_mu_sm"
                block_lbl{k} = [namePrefix '_' predName];
            end
            
            % In this version, every predictor defined is "active" 
            % so the mask is all true for the defined predictors.
            block_mask = true(1, numPreds);
            
            % Collect results
            xtrue_parts{end+1} = block_xtrue;
            lbl_parts{end+1} = block_lbl;
            mask_parts{end+1} = block_mask;
        end
    end
    
    % Concatenate all parts into flat arrays
    xtrue = horzcat(xtrue_parts{:});
    lbl = horzcat(lbl_parts{:});
    mask = horzcat(mask_parts{:});
end