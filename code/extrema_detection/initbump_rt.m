
function s = initbump_rt(varargin)


%% Parse inputs
p = inputParser;
addParameter(p, 'export', false);
addParameter(p, 'template_length', '');
addParameter(p, 'frames_2b_exp', '');
addParameter(p, 'run2analyse', []);  % 'offset' (old behavior) or 'onset'
addParameter(p, 'which_2_test', '');

parse(p, varargin{:});

export = p.Results.export;
d = p.Results.template_length;
frames_2b_exp = p.Results.frames_2b_exp;
run2analyse = p.Results.run2analyse;
test = p.Results.which_2_test;

% Specify Paths
run = sprintf('run%02d', run2analyse);
paths = path_generator('folder', fullfile('/extrema_detection', test, run));

% Load the Motion Cache
if strcmp(test, 'ground_truth') || strcmp(test, 'systematic_analysis')
    sim_params = importdata(fullfile(paths.results, 'sim_params.mat'));

elseif strcmp(test, 'empirical')
    sim_params = importdata(fullfile(paths.results, 'fit_results.mat'));

end
motion_cache = importdata(sim_params.motion_cache_path);

% Calculations
extra_frames = d - frames_2b_exp;
d_s = d/60;

% d is the template length
% Choose a weighting profile (pick one)

% 1) Linear weights: small -> large across the template
% w = linspace(0.5, 1.0, d).';              % adjust endpoints as desired

% 2) Exponential weights (stronger emphasis on latest samples)
% lambda controls how fast weights grow; try 2–5
% idx 0..d-1 so the last entries get the largest weight

% Load relevant files
file_list = dir(paths.results);
folders = file_list([file_list.isdir]);
folders = folders(~ismember({folders.name}, {'.', '..'}));

for idx_file = 1:size(folders, 1)
    fprintf('folder %d out of %d \n',idx_file, size(folders, 1))
    tic
    paths.results = fullfile(folders(idx_file).folder, folders(idx_file).name);
    y.ed = importdata(fullfile(paths.results, 'sims_ed/y.mat'));
    y.ac = importdata(fullfile(paths.results, 'sims_ac/y.mat'));

    for exp_gradient = 0
        w = exp(linspace(0, exp_gradient, d)).';            % example with growth factor ~e^2

        w = w * (d / sum(w));

        code4segm = sprintf('d%d_2bexp%d_expkern%d', d, frames_2b_exp, exp_gradient);

        % Select only freezes with specific durations
        for idx_gen_model = {'ac', 'ed'}

            gen_model = idx_gen_model{1};
            gen_path = sprintf('sims_%s', gen_model);

            paths_out = path_generator('folder', fullfile('/extrema_detection', test, run, folders(idx_file).name, gen_path, code4segm));
            mkdir(paths_out.fig); mkdir(paths_out.results);

            y_curr = y.(gen_model);
            y_curr = y_curr(y_curr.durations_s > d_s & y_curr.durations_s < sim_params.T, :);

            flies = y_curr.fly;
            allfr_onsets = y_curr.onsets;
            allfr_durations = round(y_curr.durations_s .* 60);
            allfr_offsets = allfr_onsets + allfr_durations - 1;

            s = struct;

            comparison_sm_cropped_all = cell(1, height(y_curr));

            for idx_comparison = 1:height(y_curr)
                comparison_sm = motion_cache(flies(idx_comparison));
                comparison_onset = allfr_onsets(idx_comparison);
                comparison_offset_wextra = allfr_offsets(idx_comparison) + extra_frames - 1;
                comparison_sm_cropped_all{idx_comparison} = comparison_sm(comparison_onset:comparison_offset_wextra);
            end

            parfor idx_bout = 1:height(y_curr)

                sm = motion_cache(flies(idx_bout));
                freeze_onset = allfr_onsets(idx_bout);
                freeze_offset = allfr_offsets(idx_bout);

                sm_bout = sm(freeze_onset:freeze_offset);

                template_onset = length(sm_bout) - d + 1;

                v1 = sm_bout(template_onset:end);

                bout = table();
                similarity_value = nan(1, height(y_curr));
                similarity_sort = nan(1, height(y_curr));
                summed_motion_b4 = nan(1, height(y_curr));
                summed_motion_during = nan(1, height(y_curr));
                rt_post_template = nan(1, height(y_curr));

                for idx_comparison = 1:height(y_curr)

                    comparison_sm_cropped = comparison_sm_cropped_all{idx_comparison};

                    % Precompute weighted template terms
                    v1_wsq = sum(w .* (v1.^2));                % scalar
                    dots_w = conv(comparison_sm_cropped, flipud(w .* v1), 'valid');
                    ssq_w  = conv(comparison_sm_cropped.^2, flipud(w), 'valid');

                    % Weighted SSD distance
                    similarity_framebframe_vectorized = v1_wsq + ssq_w - 2*dots_w;

                    % The rest stays the same
                    [closest_similarity, best_frame] = min(similarity_framebframe_vectorized);
                    similarity_value(idx_comparison) = closest_similarity;
                    similarity_sort(idx_comparison) = best_frame;
                    summed_motion_b4(idx_comparison) = sum(comparison_sm_cropped(1:best_frame - 1));
                    summed_motion_during(idx_comparison) = sum(comparison_sm_cropped(best_frame:best_frame + d - 1));

                    rt_post_template(idx_comparison) = allfr_durations(idx_comparison) - best_frame + 1;

                end

                bout.closest_similarity = similarity_value';
                bout.best_frame = similarity_sort';
                bout.summed_motion_b4 = summed_motion_b4';
                bout.summed_motion_during = summed_motion_during';

                bout.rt_post_template = rt_post_template';
                bout.idx_freeze = (1:height(y_curr))';

                sorted_bout = sortrows(bout, 'closest_similarity');
                s(idx_bout).boutlist = sorted_bout;
                s(idx_bout).sm = sm_bout;

                
            end
            
            size(s, 2) == height(y_curr)
            cd(paths_out.results)
            save(sprintf('struct_%s.mat', gen_model), 's')
            save(sprintf('comparison_sm_cropped_%s.mat', gen_model), 'comparison_sm_cropped_all')
        end
    end
toc
end
end