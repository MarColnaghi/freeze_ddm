function bouts = generate_bouts_table(gap_imm, varargin)
% GENERATE_BOUTS_TABLE  Builds the bouts table for one fill-in setting and
% saves it where path_generator looks for it.
%
%   bouts = generate_bouts_table(3)
%   bouts = generate_bouts_table(2, 'imfirst', true)
%   bouts = generate_bouts_table(3, 'dataset', 'LC11')
%   for gap = 0:3, generate_bouts_table(gap); end        % a sweep
%
% gap_imm is the number in the folder name: imm3_mob3_pc4 bridges moving gaps
% of up to 3 frames and drops immobile runs of up to 3 frames. In the
% thresholds that is fill_in = gap + 1, because bwareaopen(x, P) removes runs
% shorter than P frames. The table itself is built in load_flies_new.
%
% Options
%   gap_mob   gap_imm    longest moving gap bridged into a bout (imm and mob
%                        were swept together in the old driver, hence the default)
%   pc        4          pixelchange threshold, the pc in the folder name
%   imfirst   false      true removes small immobile runs before filling
%                        moving gaps, and uses fill_in_combinations_imfirst
%   overwrite false      an existing bouts.mat is kept, and loaded if an output
%                        is asked for. true rebuilds it. The old file is only
%                        replaced once the new table is complete
%   save      'bouts'    'both' also saves soc_mot.mat, the motion trace of
%                        every bout (large)
%   out_dir   ''         save here instead of the folder from path_generator,
%                        to try a setting without touching datasets/
%   dataset   'CS'       which experiment to read, see load_flies_new: 'CS'
%                        (004--social_ddm) or 'LC11' (006--lc11). LC11 saves
%                        under datasets/LC11/, so its tables do not collide
%                        with the CS ones of the same id
%
% A folder this call created is removed again if the build fails before
% anything was written to it, so a failed run leaves no empty setting behind.
%
% The freeze_bout label in dataset_resolved, and so freeze_cache.mat, is built
% at gap 3 (params.freeze_gap in generate_dataset_resolved) to match imm3. For
% any other setting it will not match this table unless freeze_gap follows.
% That is CS only, LC11 has no dataset_resolved or caches built yet.

opt = inputParser;
addRequired(opt, 'gap_imm', @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}));
addParameter(opt, 'gap_mob', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0 && x == round(x)));
addParameter(opt, 'pc', 4, @(x) isnumeric(x) && isscalar(x));
addParameter(opt, 'imfirst', false, @(x) islogical(x) && isscalar(x));
addParameter(opt, 'overwrite', false, @(x) islogical(x) && isscalar(x));
addParameter(opt, 'save', 'bouts', @(x) any(strcmp(x, {'bouts', 'both'})));
addParameter(opt, 'out_dir', '', @(x) ischar(x) || isstring(x));
addParameter(opt, 'dataset', 'CS', @(x) any(strcmp(x, {'CS', 'LC11'})));
parse(opt, gap_imm, varargin{:});

gap_mob = opt.Results.gap_mob;
if isempty(gap_mob)
    gap_mob = gap_imm;
end

thresholds = define_thresholds();
thresholds.fill_in_imm = gap_imm + 1;
thresholds.fill_in_mob = gap_mob + 1;
thresholds.pc = opt.Results.pc;

id_code = sprintf('imm%d_mob%d_pc%d', gap_imm, gap_mob, opt.Results.pc);
dataset_folder = opt.Results.dataset;
if strcmp(dataset_folder, 'CS')
    dataset_folder = '';
end
paths = path_generator('bouts_id', id_code, 'imfirst', opt.Results.imfirst, 'dataset', dataset_folder);
if strlength(opt.Results.out_dir) > 0
    paths.dataset = char(opt.Results.out_dir);
end
bouts_file = fullfile(paths.dataset, 'bouts.mat');

if isfile(bouts_file) && ~opt.Results.overwrite
    fprintf('already have %s in %s, not rebuilding (pass ''overwrite'', true to redo it)\n', id_code, paths.dataset)
    if nargout > 0
        bouts = importdata(bouts_file);
    end
    return
end

made_dir = ~isfolder(paths.dataset);
if made_dir
    mkdir(paths.dataset);
end

tic
try
    bouts = load_flies_new(thresholds, 'paths', paths, 'save', opt.Results.save, 'imfirst', opt.Results.imfirst, 'dataset', opt.Results.dataset);
catch err
    if made_dir && isempty(setdiff({dir(paths.dataset).name}, {'.', '..'}))
        rmdir(paths.dataset);
    end
    rethrow(err)
end
fprintf('built %s in %.1f min, saved to %s\n', id_code, toc / 60, paths.dataset)
end
