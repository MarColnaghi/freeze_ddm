function eval_fun = compile_evaluate_model(model, y, lbl)
% COMPILE_EVALUATE_MODEL
% Builds a fast evaluator for evaluate_model
%
% model : original model struct
% y     : STRUCT with predictor fields (not table!)
% lbl   : parameter labels from get_ground_truth_vector

    param_names = fieldnames(model);
    blocks = struct([]);
    beta_idx = 1;

    for i = 1:numel(param_names)
        pname = param_names{i};
        block = model.(pname);

        blocks(i).name = pname;
        blocks(i).link = block.link;

        preds = block.predictors;
        nPred = numel(preds);

        blocks(i).X = cell(nPred,1);
        blocks(i).is_cell = false(nPred,1);
        blocks(i).beta_inds = zeros(nPred,1);

        for j = 1:nPred
            pred_name = preds{j}.name;

            % Store predictor data directly
            if ~isfield(y, pred_name)
                error('Predictor "%s" not found.', pred_name);
            end

            blocks(i).X{j} = y.(pred_name);
            blocks(i).is_cell(j) = iscell(y.(pred_name));

            % Find corresponding beta index ONCE
            tf = strcmp(lbl, [pname '_' pred_name]);
            if ~any(tf)
                error('Coefficient %s_%s not found.', pname, pred_name);
            end

            blocks(i).beta_inds(j) = find(tf,1);
        end
    end

    % Return fast evaluator
    eval_fun = @(x) eval_compiled_model(blocks, x);
end
