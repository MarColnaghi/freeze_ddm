function check_model_consistency(model)
% CHECK_MODEL_CONSISTENCY validates that all fields of the model are consistent.
%
% Throws an error if inconsistencies are found.
% Checks include: fn validity, param index size, bounds format, etc.

fields = fieldnames(model);
num_params_used = [];

for i = 1:numel(fields)
    key = fields{i};
    spec = model.(key);

    % Check required fields exist
    required = {'params', 'fn', 'bounds'};
    for r = required
        if ~isfield(spec, r{1})
            error('Field "%s" is missing required subfield "%s".', key, r{1});
        end
    end

    % Check fn is function handle
    if ~isa(spec.fn, 'function_handle')
        error('Field "%s": fn is not a function handle.', key);
    end

    % Check bounds format
    b = spec.bounds;
    if ~isnumeric(b) || ~isequal(size(b), [1, 4])
        error('Field "%s": bounds must be a 1x4 numeric array.', key);
    end

    % Check params are numeric indices
    p_indices = spec.params;
    if ~isnumeric(p_indices)
        error('Field "%s": params must be numeric indices.', key);
    end

    % Track all indices for duplicate detection later
    num_params_used = [num_params_used, p_indices];

    % Try evaluating the function with dummy input
    try
        dummy_p = ones(1, numel(p_indices));  % 1 for each param
        dummy_y = struct();  % you can expand this later if needed
        out = spec.fn(dummy_p, dummy_y);
        
        % Optionally: check shape or type of output
        if ~isnumeric(out)
            warning('Field "%s": fn output is not numeric.', key);
        end
    catch ME
        error('Field "%s": error when evaluating fn: %s', key, ME.message);
    end
end

% Optional: check for duplicate param indices
if numel(unique(num_params_used)) ~= numel(num_params_used)
    warning('Some parameter indices are reused across fields.');
end

disp('Model passed all consistency checks.');

end
