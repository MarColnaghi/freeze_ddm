function [xtrue, lbl, mask] = get_ground_truth_vector(model)
    
components = fieldnames(model);
    xtrue_parts = cell(1, numel(components));
    xtrue_labels = cell(1, numel(components));
    xtrue_mask = cell(1, numel(components));
    
    base_names = {'sm', 'smp', 'fs', 'ln', 'ls', 'tsl', 'intercept'};
    
    for i = 1:numel(components)
        comp = model.(components{i});
        vec = nan(1, numel(base_names));
        
        % Check if the ground_truth field exists
        has_gt = isfield(comp, 'ground_truth');
        
        for idx_predictors = 1:size(comp.predictors, 2)
            % Identify which base_name this predictor matches
            selector = strcmp(comp.predictors{idx_predictors}.name, base_names);
            
            if any(selector)
                if has_gt
                    % Use actual ground truth if it exists
                    vec(selector) = comp.ground_truth(idx_predictors);
                else
                    % Default to zero as requested
                    vec(selector) = 0;
                end
            end
        end
        
        xtrue_parts{i} = vec;
        xtrue_labels{i} = strcat(components{i}, "_", base_names);
        xtrue_mask{i} = ~isnan(vec);
    end
    
    xtrue = horzcat(xtrue_parts{:});
    lbl = horzcat(xtrue_labels{:});
    mask = horzcat(xtrue_mask{:});
end