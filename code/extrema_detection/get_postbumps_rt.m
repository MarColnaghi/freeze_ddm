function [distance, peaks] = get_postbumps_rt(varargin)


%% Parse inputs
p = inputParser;
addParameter(p, 'export', false);
addParameter(p, 'template_length', '');
addParameter(p, 'frames_2b_exp', '');
addParameter(p, 'n_selected_comparisons', '');
addParameter(p, 'gen_model', '');

addParameter(p, 'run2analyse', []); 

parse(p, varargin{:});

export = p.Results.export;
d = p.Results.template_length;
frames_2b_exp = p.Results.frames_2b_exp;
n_selected_comparisons = p.Results.n_selected_comparisons;
run2analyse = p.Results.run2analyse;
gen_model = p.Results.gen_model;

% Specify Paths
run = sprintf('run%02d', run2analyse);
paths = path_generator('folder', fullfile('/extrema_detection', 'ac_vs_ed', run));

% Load the Motion Cache
sim_params = importdata(fullfile(paths.results, 'sim_params.mat'));
motion_cache = importdata(sim_params.motion_cache_path);

% Calculations
extra_frames = d - frames_2b_exp - 1;
d_s = d/60;

% Select Color Map
col = cbrewer2('Spectral', n_selected_comparisons);

% Already Creates the figure for the histogram
fh = figure('color', 'w', 'Position', [100, 100, 600, 250]);
hold on

% Load relevant files
y.ed = importdata(fullfile(paths.results, 'sims_ed/y.mat'));
y.ac = importdata(fullfile(paths.results, 'sims_ac/y.mat'));

code4segm = sprintf('d_%d-2bexp_%d-ncomp_%d', d, frames_2b_exp, n_selected_comparisons);
paths_out = path_generator('folder', fullfile('/extrema_detection', 'ac_vs_ed', run, code4segm));
mkdir(paths_out.fig); mkdir(paths_out.results);

% Select only freezes with specific durations
y_curr = y.(gen_model);
y_curr = y_curr(y_curr.durations_s < sim_params.T, :);

flies = y_curr.fly;
allfr_onsets = y_curr.onsets;
allfr_durations = round(y_curr.durations_s .* 60) + 1;
allfr_offsets = allfr_onsets + allfr_durations;

s = struct;

fh = figure('color','w','Position',[100, 100, 1200, 650]);
ax = gca;
apply_generic(ax, 18)
xlabel('Time aligned to Onset (frames)')
ax.YAxis.Visible = 'off';
hold on
distance = nan(1, height(y_curr));
peaks = nan(1, height(y_curr));

for idx_bout = 1:height(y_curr)

    fprintf('bout n:%d\n', idx_bout);
    sm = motion_cache(flies(idx_bout));
    
    freeze_onset = allfr_onsets(idx_bout);
    freeze_duration = allfr_durations(idx_bout);
    freeze_offset = allfr_offsets(idx_bout);

    sm_bout = sm(freeze_onset:freeze_offset);

    [pks, idx] = findpeaks(sm_bout, 'MinPeakHeight', 1);

    if idx_bout < 80
        plh = plot(sm_bout + idx_bout, 'k');
        scatter(idx, [pks + idx_bout + 0.3], 'MarkerFaceColor', 'r', 'Marker', 'v', 'MarkerEdgeColor', 'none', 'SizeData', 8)
    end
    if ~isempty(pks)
        peaks(idx_bout) = pks(end);
        distance(idx_bout) = idx(end) - freeze_duration;
    end

end

figure
histogram(distance, -600.5:5:1.5)
ylabel('Count')
xlabel('Distance from Last Peak (frames)')

apply_generic(gca)
figure
scatter(distance, y_curr.sm)