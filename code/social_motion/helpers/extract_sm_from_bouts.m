
function sm_output = extract_sm_from_bouts(bouts_le, varargin)

opt = inputParser;
addParameter(opt, 'type', 'onsets');
addParameter(opt, 'window', [-180 300]);
addParameter(opt, 'output_type', 'mat');
addParameter(opt, 'size', 630);
addParameter(opt, 'align', 'onset');
addParameter(opt, 'delay_after', 0);
addParameter(opt, 'cache', 'motion_cache');
addParameter(opt, 'norm_factor', 1);

parse(opt, varargin{:});

norm_factor = opt.Results.norm_factor;
delay_after = opt.Results.delay_after;
cache = opt.Results.cache;
window = opt.Results.window;

paths = path_generator('folder', 'social_motion');
cache = importdata(fullfile(paths.cache_path, strcat(cache, '.mat')));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));


switch opt.Results.type

    case 'loom'

        offsets = (window(1) : window(2));
        total_looms = 20;
        total_flies = size(cache, 1);

        n_moving_flies = nan(total_flies, 1);
        sloom = nan(total_flies, 1);
        sm_output = nan(total_flies, total_looms, length(offsets));

        for idx_fly = 1:total_flies

            sm_fly = cache(idx_fly);
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
            sum_motion = cache(fly_idx);

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
                    sm_output{idx_bouts} = chunk_data ./ norm_factor;

                case 'mat'
                    sm_output(idx_bouts, :) = chunk_data ./ norm_factor;
            end
        end

    case 'onlyfreeze'
        total_bouts = height(bouts_le);
        switch opt.Results.output_type
            case 'cell'
                sm_output = cell(total_bouts, 1);
                for idx_bouts = 1:total_bouts
                    sum_motion = cache(bouts_le.fly(idx_bouts));
                    max_idx = numel(sum_motion);

                    ons = bouts_le.onsets(idx_bouts);
                    % Consistent calculation of requested end
                    req_off = bouts_le.ends(idx_bouts) - 1 + delay_after;
                    off = min(req_off, max_idx);

                    % Force column vector (:) to prevent concatenation errors
                    sig = sum_motion(ons:off) ./ norm_factor;

                    if req_off > max_idx
                        padding_needed = req_off - max_idx;
                        sig = [sig(:); nan(padding_needed, 1)];
                    end
                    sm_output{idx_bouts} = sig;
                end

            case 'mat'
                % FIX 1: Ensure max_len exactly matches the logic inside the loop
                % We calculate the width based on the intended window size for all rows
                freeze_lens = (bouts_le.ends - bouts_le.onsets) + delay_after;
                max_len = max(freeze_lens);
                sm_output = nan(total_bouts, max_len);

                for idx_bouts = 1:total_bouts
                    sum_motion = cache(bouts_le.fly(idx_bouts));
                    max_idx = numel(sum_motion);

                    ons = bouts_le.onsets(idx_bouts);
                    req_off = bouts_le.ends(idx_bouts) - 1 + delay_after;
                    off = min(req_off, max_idx);

                    % Extract the available signal
                    freeze_sig = sum_motion(ons:off) ./ norm_factor;

                    % FIX 2: Create a row vector of the exact width needed for this bout
                    % This ensures it matches the pre-calculated 'freeze_lens'
                    this_bout_len = freeze_lens(idx_bouts);
                    full_segment = nan(1, this_bout_len);

                    % Place the actual signal into the segment (the rest remains NaN)
                    L_actual = numel(freeze_sig);
                    full_segment(1:L_actual) = freeze_sig;

                    % FIX 3: Robust Assignment
                    % We use 1:this_bout_len to ensure we don't exceed max_len
                    switch opt.Results.align
                        case 'onset'
                            sm_output(idx_bouts, 1:this_bout_len) = full_segment;
                        case 'offset'
                            % Right-aligning within the matrix
                            start_col = max_len - this_bout_len + 1;
                            sm_output(idx_bouts, start_col:end) = full_segment;
                    end
                end
        end


end


