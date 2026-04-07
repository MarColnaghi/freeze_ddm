
function sm_output = extract_sm_from_bouts(bouts_le, varargin)

opt = inputParser;
addParameter(opt, 'type', 'onsets');
addParameter(opt, 'window', [-180 300]);
addParameter(opt, 'output_type', 'mat');
addParameter(opt, 'size', 630);
addParameter(opt, 'align', 'onset');
addParameter(opt, 'delay', 0);

parse(opt, varargin{:});

delay = opt.Results.delay;

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

        win_start = opt.Results.window(1);
        win_end   = opt.Results.window(2);
        win_width = opt.Results.window(2) - opt.Results.window(1) + 1;
        
        switch opt.Results.output_type
            case 'cell'
                % Pre-allocate a cell array for speed
                sm_output = cell(total_bouts, 1);

            case 'mat'
                % Pre-allocate a matrix filled with NaNs (or zeros)
                % This is significantly faster than growing the matrix inside the loop
                sm_output = nan(total_bouts, win_width);
        end

        for idx_bouts = 1:total_bouts
            ons = bouts_le.onsets(idx_bouts);
            fly_idx = bouts_le.fly(idx_bouts);
            sum_motion = motion_cache(fly_idx);

            % Calculate absolute indices in the motion data
            idx_range = (ons + win_start) : (ons + win_end);

            % Boundary Check: Ensure indices are within the actual data length
            % This prevents "Index out of bounds" errors for early or late onsets
            valid_idx = idx_range > 0 & idx_range <= length(sum_motion);

            % Initialize a NaN or zero vector of the required window size
            chunk_data = nan(1, length(idx_range));

            % Fill only the valid portions
            chunk_data(valid_idx) = sum_motion(idx_range(valid_idx));

            % Normalize and assign
            switch opt.Results.output_type
                case 'cell'
                    sm_output{idx_bouts} = chunk_data ./ 10;

                case 'mat'
                    sm_output(idx_bouts, :) = chunk_data ./ 10;
            end
        end

    case 'onlyfreeze'

        total_bouts = height(bouts_le);

        switch opt.Results.output_type

            case 'cell'

                sm_output = cell(total_bouts, 1);

                for idx_bouts = 1:total_bouts
                    ons = bouts_le.onsets(idx_bouts) + delay;
                    off = bouts_le.ends(idx_bouts) - 1;
                    sum_motion = motion_cache(bouts_le.fly(idx_bouts));

                    sm_output{idx_bouts} = sum_motion(ons:off) ./ 10;
                end

            case 'mat'

                freeze_lens = bouts_le.ends - bouts_le.onsets;
                max_len = max(freeze_lens);

                sm_output = nan(total_bouts, max_len);

                for idx_bouts = 1:total_bouts
                    ons = bouts_le.onsets(idx_bouts) + delay;
                    off = bouts_le.ends(idx_bouts) - 1;
                    
                    sum_motion = motion_cache(bouts_le.fly(idx_bouts));

                    freeze_sig = sum_motion(ons:off) ./ 10;
                    L = numel(freeze_sig);

                    switch opt.Results.align
                        case 'onset'
                            
                            sm_output(idx_bouts, 1:L) = freeze_sig;
                        case 'offset'
                            sm_output(idx_bouts, end-L+1:end) = freeze_sig;
                    end

                    
                    
                end

        end


end


