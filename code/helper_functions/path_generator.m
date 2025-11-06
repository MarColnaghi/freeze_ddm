function [path] = path_generator(varargin)

opt = inputParser;
addParameter(opt, 'bouts_id', 'imm3_mob3_pc4');
addParameter(opt, 'folder', []);
addParameter(opt, 'model', []);

parse(opt, varargin{:});

% Define paths
if nargin > 1
    path.results = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results', opt.Results.folder, opt.Results.model);
    path.filename.elbo = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/', opt.Results.folder, opt.Results.model, 'elbo.mat');
    path.filename.est_params = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results', opt.Results.folder, opt.Results.model, 'est_params.mat');
    path.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures', opt.Results.folder, opt.Results.model);
else
    path.results = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results', opt.Results.folder);
    path.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures', opt.Results.folder);
end

path.code = '/Users/marcocolnaghi/PhD/freeze_ddm/code';
path.dataset = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/datasets/fill_in_combinations_mobfirst/', opt.Results.bouts_id);
path.cache_path = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/datasets/caches');

