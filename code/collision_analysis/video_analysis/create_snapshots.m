% create_snapshots.m
% Superposition analysis: builds per-session min projections and global
% mean / min (occupancy) images at several inter-fly distance thresholds.
%
% Requires: Image Processing Toolbox (rgb2gray, imresize, imwrite)

% --- CONFIGURATION ---
THRESHOLDS           = [40, 50, 55, 60, 65, 70, 80, 100];
EPSILON              = 1.0;
CROP_SIZE            = 200;
SAMPLES_PER_PROJ     = 5;
VIDEO_DIR_BASE       = '/Volumes/Elements/SocialDDM-CS/aligned';
DIST_DATA_DIR        = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_0';
OUTPUT_DIR           = 'superposition_results';

% ---------------------

run_superposition_organized(THRESHOLDS, EPSILON, CROP_SIZE, SAMPLES_PER_PROJ, ...
    VIDEO_DIR_BASE, DIST_DATA_DIR, OUTPUT_DIR);

% =========================================================

function run_superposition_organized(THRESHOLDS, EPSILON, CROP_SIZE, SAMPLES_PER_PROJ, ...
        VIDEO_DIR_BASE, DIST_DATA_DIR, OUTPUT_DIR)

n_thresh  = length(THRESHOLDS);
half_crop = CROP_SIZE / 2;

% --- Setup directories ---
proj_base_dir = fullfile(OUTPUT_DIR, 'min_projections');
if ~exist(proj_base_dir, 'dir'); mkdir(proj_base_dir); end
for i = 1:n_thresh
    d = fullfile(proj_base_dir, sprintf('dist_%d', THRESHOLDS(i)));
    if ~exist(d, 'dir'); mkdir(d); end
end

% --- Find distance CSV files ---
dist_listing = dir(fullfile(DIST_DATA_DIR, '*.fly*.csv'));
fprintf('Found %d distance files.\n', length(dist_listing));

% --- Global accumulators ---
accum_sum  = cell(1, n_thresh);
accum_min  = cell(1, n_thresh);
counters   = zeros(1, n_thresh);
for i = 1:n_thresh
    accum_sum{i} = zeros(CROP_SIZE, CROP_SIZE, 'single');
    accum_min{i} = uint8(255 * ones(CROP_SIZE, CROP_SIZE));
end

for df = 1:length(dist_listing)
    dist_path = fullfile(dist_listing(df).folder, dist_listing(df).name);
    filename  = dist_listing(df).name;

    % Extract fly_id from filename (e.g. "session.fly3.csv" -> "fly3")
    parts  = strsplit(filename, '.');
    fly_id = '';
    for p = 1:length(parts)
        if contains(parts{p}, 'fly')
            fly_id = parts{p};
            break;
        end
    end
    if isempty(fly_id); continue; end

    fly_folder   = fullfile(VIDEO_DIR_BASE, fly_id);
    coord_listing = dir(fullfile(fly_folder, 'aligned_*_coordinates.csv'));

    for cf = 1:length(coord_listing)
        coord_path = fullfile(coord_listing(cf).folder, coord_listing(cf).name);
        session_id = strrep(strrep(coord_listing(cf).name, 'aligned_', ''), '_coordinates.csv', '');
        video_path = fullfile(fly_folder, sprintf('aligned_%s.mp4', session_id));

        if ~exist(video_path, 'file'); continue; end

        % Skip session if all threshold images already exist
        needed = false(1, n_thresh);
        for i = 1:n_thresh
            out_path = fullfile(proj_base_dir, sprintf('dist_%d', THRESHOLDS(i)), ...
                sprintf('%s_%s_minproj.png', fly_id, session_id));
            if ~exist(out_path, 'file')
                needed(i) = true;
            end
        end
        if ~any(needed)
            fprintf('Skipping %s: all projection images already exist.\n', session_id);
            continue;
        end

        try
            dist_df    = readtable(dist_path);
            coord_df   = readtable(coord_path, 'ReadVariableNames', false);
            distances  = dist_df.dist_min;
            % Coordinates are stored as 0-indexed pixel values (OpenCV convention)
            focal_coords = coord_df{:, 1:2};
            min_len    = min(length(distances), size(focal_coords, 1));

            % Find 1-indexed frame positions for each threshold
            thresh_idx = cell(1, n_thresh);
            for i = 1:n_thresh
                thresh_idx{i} = find(abs(distances(1:min_len) - THRESHOLDS(i)) < EPSILON);
            end

            % Random subset for per-session projections
            frames_to_sample = cell(1, n_thresh);
            all_needed_frames = [];
            for i = 1:n_thresh
                avail = thresh_idx{i};
                if ~isempty(avail)
                    k = min(length(avail), SAMPLES_PER_PROJ);
                    sel = avail(randperm(length(avail), k));
                    frames_to_sample{i} = sel;
                    all_needed_frames   = union(all_needed_frames, avail);
                end
            end
            if isempty(all_needed_frames); continue; end

            needed_set = false(min_len, 1);
            needed_set(all_needed_frames(all_needed_frames <= min_len)) = true;

            vid      = VideoReader(video_path);
            max_frame = max(all_needed_frames);
            curr_idx = 0;
            proj_bufs = cell(1, n_thresh);

            fprintf('Processing %s (streaming to frame %d)...\n', session_id, max_frame);

            while hasFrame(vid) && curr_idx < max_frame
                frame    = readFrame(vid);
                curr_idx = curr_idx + 1;

                if curr_idx <= min_len && needed_set(curr_idx)
                    % Convert 0-indexed OpenCV coords to MATLAB 1-indexed slice bounds
                    cx  = focal_coords(curr_idx, 1);
                    cy  = focal_coords(curr_idx, 2);
                    y_s = cy - half_crop + 1;
                    y_e = cy + half_crop;
                    x_s = cx - half_crop + 1;
                    x_e = cx + half_crop;

                    if y_s >= 1 && x_s >= 1 && y_e <= size(frame,1) && x_e <= size(frame,2)
                        gray_crop = rgb2gray(frame(y_s:y_e, x_s:x_e, :));

                        for i = 1:n_thresh
                            if ismember(curr_idx, thresh_idx{i})
                                accum_sum{i}  = accum_sum{i} + single(gray_crop);
                                accum_min{i}  = min(accum_min{i}, gray_crop);
                                counters(i)   = counters(i) + 1;
                            end
                            if ismember(curr_idx, frames_to_sample{i})
                                proj_bufs{i}{end+1} = gray_crop;
                            end
                        end
                    end
                end
            end
            clear vid;

            % Save per-session min projections
            for i = 1:n_thresh
                if ~isempty(proj_bufs{i})
                    stack    = cat(3, proj_bufs{i}{:});
                    min_proj = min(stack, [], 3);
                    out_path = fullfile(proj_base_dir, sprintf('dist_%d', THRESHOLDS(i)), ...
                        sprintf('%s_%s_minproj.png', fly_id, session_id));
                    imwrite(normalize_minmax(min_proj), out_path);
                end
            end

        catch ME
            fprintf('Error processing %s: %s\n', session_id, ME.message);
        end
    end
end

% --- Save global results ---
fprintf('\nSaving global aggregated results...\n');
for i = 1:n_thresh
    if counters(i) > 0
        global_mean_norm = normalize_minmax(uint8(accum_sum{i} / counters(i)));
        imwrite(global_mean_norm, fullfile(OUTPUT_DIR, sprintf('global_mean_dist_%d.png', THRESHOLDS(i))));

        global_min_norm  = normalize_minmax(accum_min{i});
        imwrite(global_min_norm,  fullfile(OUTPUT_DIR, sprintf('global_min_occupancy_dist_%d.png', THRESHOLDS(i))));
    end
end
fprintf('Task completed!\n');
end

% ---------------------------------------------------------

function out = normalize_minmax(img)
    img_d = double(img);
    mn = min(img_d(:));
    mx = max(img_d(:));
    if mx > mn
        out = uint8(255 * (img_d - mn) / (mx - mn));
    else
        out = uint8(img_d);
    end
end
