function [bouts_sloom, mask] = quantilizer(bouts, varargin)

% Parse inputs
opt = inputParser;
addParameter(opt, 'quant', struct('sm', 3, 'fs', 2));  % Default quant values
addParameter(opt, 'idx_quanti', []);
parse(opt, varargin{:});

quant = opt.Results.quant;
idx_quanti = opt.Results.idx_quanti;

% Extract index fields
if isfield(idx_quanti, 'ln')
    [idx_loom_speeds, idx_speeds, idx_motion, idx_nloom] = ...
        deal(idx_quanti.ls, idx_quanti.fs, idx_quanti.sm, idx_quanti.ln);
else
    [idx_loom_speeds, idx_speeds, idx_motion] = ...
        deal(idx_quanti.ls, idx_quanti.fs, idx_quanti.sm);
end

% Select Loom
mask_sloom = bouts.sloom_norm == idx_loom_speeds;
bouts_sloom = bouts(mask_sloom, :);
first_mask = find(mask_sloom);

% Compute quantiles and discretize
thr_sm = quantile(bouts_sloom.avg_sm_freeze_norm, linspace(0, 1, quant.sm + 1));
quant_sm = discretize(bouts_sloom.avg_sm_freeze_norm, thr_sm);

thr_fs = quantile(bouts_sloom.avg_fs_1s_norm, linspace(0, 1, quant.fs + 1));
quant_fs = discretize(bouts_sloom.avg_fs_1s_norm, thr_fs);

% Final selection
if isfield(idx_quanti, 'ln')
    thr_ln = 0:5:20;
    quant_ln = discretize(bouts_sloom.nloom, thr_ln);
    bouts_sloom = bouts_sloom(quant_sm == idx_motion & ...
        quant_fs == idx_speeds & ...
        quant_ln == idx_nloom, :);

    mask = first_mask(quant_sm == idx_motion & ...
        quant_fs == idx_speeds & ...
        quant_ln == idx_nloom);
else
    bouts_sloom = bouts_sloom(quant_sm == idx_motion & ...
        quant_fs == idx_speeds, :);

    mask = first_mask(quant_sm == idx_motion & ...
        quant_fs == idx_speeds);
end
