function bouts = quantilizer(bouts, varargin)

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
bouts = bouts(bouts.sloom_norm == idx_loom_speeds,:);

% Compute quantiles and discretize
thr_sm = quantile(bouts.avg_sm_freeze_norm, linspace(0, 1, quant.sm + 1));
quant_sm = discretize(bouts.avg_sm_freeze_norm, thr_sm);

thr_fs = quantile(bouts.avg_fs_1s_norm, linspace(0, 1, quant.fs + 1));
quant_fs = discretize(bouts.avg_fs_1s_norm, thr_fs);

% Final selection
if isfield(idx_quanti, 'ln')
    thr_ln = 0:5:20;
    quant_ln = discretize(bouts.nloom, thr_ln);
    bouts = bouts(quant_sm == idx_motion & ...
                  quant_fs == idx_speeds & ...
                  quant_ln == idx_nloom,:);
else
    bouts = bouts(quant_sm == idx_motion & ...
                  quant_fs == idx_speeds,:);
end
