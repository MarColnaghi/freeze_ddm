
function sm_output = extract_sm_from_bouts(bouts_le, varargin)

opt = inputParser;
addParameter(opt, 'type', 'onsets');
addParameter(opt, 'window', [-180 300]);
addParameter(opt, 'formula', 'mat');
addParameter(opt, 'size', 630);
addParameter(opt, 'align', 'onset');

parse(opt, varargin{:});

chunk_len = opt.Results.size;
window = opt.Results.window;

paths = path_generator('folder', 'social_motion');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));


switch opt.Results.type

    case 'loom'

        offsets = (window(1) : window(2));
        total_looms = 20;
        total_flies = size(motion_cache, 1);

        n_moving_flies = nan(total_flies, 1);
        sloom = nan(total_flies, 1);
        sm_output = nan(total_flies, total_looms, length(offsets));

        for idx_fly = 1:total_flies
            sm_fly = motion_cache(idx_fly);
            freeze_frames = find(diff(loom_cache(idx_fly)) == 1);
            idx_slice = freeze_frames(:) + offsets;
            fly_loom_x_sm = sm_fly(idx_slice);
            n_moving_flies(idx_fly) = unique(bouts_le.moving_flies(bouts_le.fly == idx_fly));
            sloom(idx_fly) = unique(bouts_le.sloom(bouts_le.fly == idx_fly));

            sm_output(idx_fly, :, :) = fly_loom_x_sm;

        end


    case 'onsets'

        total_bouts = height(bouts_le);

        for idx_bouts = 1:total_bouts

            ons = bouts_le.onsets(idx_bouts);
            off = bouts_le.ends(idx_bouts) - 1;
            sum_motion = motion_cache(bouts_le.fly(idx_bouts));

            switch opt.Results.formula

                case 'cell'
                    sm_output{idx_bouts} = sum_motion(ons:ons + chunk_len) ./ 10;

                case 'mat'
                    sm_output(idx_bouts, :) = sum_motion(ons:ons + chunk_len) ./ 10;
            end


        end

    case 'onlyfreeze'

        total_bouts = height(bouts_le);

        switch opt.Results.formula

            case 'cell'

                sm_output = cell(total_bouts, 1);

                for idx_bouts = 1:total_bouts
                    ons = bouts_le.onsets(idx_bouts);
                    off = bouts_le.ends(idx_bouts) - 1;
                    sum_motion = motion_cache(bouts_le.fly(idx_bouts));

                    sm_output{idx_bouts} = sum_motion(ons:off) ./ 10;
                end

            case 'mat'

                freeze_lens = bouts_le.ends - bouts_le.onsets;
                max_len = max(freeze_lens);

                sm_output = nan(total_bouts, max_len);

                for idx_bouts = 1:total_bouts
                    ons = bouts_le.onsets(idx_bouts);
                    off = bouts_le.ends(idx_bouts) - 1;
                    sum_motion = motion_cache(bouts_le.fly(idx_bouts));

                    freeze_sig = sum_motion(ons:off) ./ 10;
                    sm_output(idx_bouts, 1:numel(freeze_sig)) = freeze_sig;
                end
        end


end


