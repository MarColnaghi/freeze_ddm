% create_snapshots_at_break.m
% Builds a global min projection from frames at the END of each freeze bout
% (the "break" frame), using a bout table as the source of frame indices.
%
% For each break frame the fly is padded, background-anchored to 255, and
% accumulated into a running minimum projection.
%
% Requires: Image Processing Toolbox (rgb2gray, imwrite)
%           Statistics and Machine Learning Toolbox (prctile)

% --- CONFIGURATION ---
CROP_SIZE       = 200;
VIDEO_DIR_BASE  = '/Volumes/Elements/SocialDDM-CS/aligned';
BOUTS_DATA_PATH = 'bouts/bouts_no_contacts_50.csv';
OUTPUT_DIR      = 'superposition_results_breaks/no_contacts_50';
% ---------------------

run_break_superposition(CROP_SIZE, VIDEO_DIR_BASE, BOUTS_DATA_PATH, OUTPUT_DIR);

% =========================================================

function run_break_superposition(CROP_SIZE, VIDEO_DIR_BASE, BOUTS_DATA_PATH, OUTPUT_DIR)

if ~exist(OUTPUT_DIR, 'dir'); mkdir(OUTPUT_DIR); end

if ~exist(BOUTS_DATA_PATH, 'file')
    fprintf('Error: %s not found.\n', BOUTS_DATA_PATH);
    return;
end

bouts_df  = readtable(BOUTS_DATA_PATH);
accum_min = uint8(255 * ones(CROP_SIZE, CROP_SIZE));
counter   = 0;
half_crop = CROP_SIZE / 2;

fly_ids = unique(bouts_df.fly);

for fi = 1:length(fly_ids)
    fly_id     = fly_ids(fi);
    fly_folder = fullfile(VIDEO_DIR_BASE, sprintf('fly%d', fly_id));
    if ~exist(fly_folder, 'dir'); continue; end

    all_files  = dir(fly_folder);
    all_names  = {all_files.name};

    video_files = all_names(cellfun(@(x) endsWith(x, '.mp4') && contains(x, 'aligned'), all_names));
    coord_files = all_names(cellfun(@(x) endsWith(x, '_coordinates.csv'), all_names));

    if isempty(video_files) || isempty(coord_files); continue; end

    video_path = fullfile(fly_folder, video_files{1});
    coord_path = fullfile(fly_folder, coord_files{1});

    try
        coord_df     = readtable(coord_path, 'ReadVariableNames', false);
        focal_coords = coord_df{:, 1:2};  % 0-indexed OpenCV pixel coords

        group        = bouts_df(bouts_df.fly == fly_id, :);
        % Python stored frame indices as 0-based; subtract 1 to get last frame
        break_frames = group.ends - 1;

        vid = VideoReader(video_path);
        fps = vid.FrameRate;

        for bi = 1:length(break_frames)
            f_idx   = break_frames(bi);   % 0-based frame index (from Python)
            mat_idx = f_idx + 1;          % MATLAB 1-based row in focal_coords

            if mat_idx > size(focal_coords, 1); continue; end

            % Seek to this frame
            vid.CurrentTime = f_idx / fps;
            if ~hasFrame(vid); continue; end
            frame = readFrame(vid);

            gray_frame = rgb2gray(frame);
            cx = focal_coords(mat_idx, 1);
            cy = focal_coords(mat_idx, 2);

            % --- Dynamic padding with median background ---
            H         = size(gray_frame, 1);
            W         = size(gray_frame, 2);
            frame_bg  = round(median(double(gray_frame(:))));
            padded    = uint8(frame_bg * ones(H + 2*half_crop, W + 2*half_crop));
            padded(half_crop+1:half_crop+H, half_crop+1:half_crop+W) = gray_frame;

            % Crop centred on fly (coordinates shifted into padded space)
            nx  = cx + half_crop;
            ny  = cy + half_crop;
            y_s = ny - half_crop + 1;
            y_e = ny + half_crop;
            x_s = nx - half_crop + 1;
            x_e = nx + half_crop;
            crop = padded(y_s:y_e, x_s:x_e);

            % --- Background anchoring: scale so 95th percentile -> 255 ---
            bg_level = prctile(double(crop(:)), 95);
            if bg_level > 0
                norm_crop = uint8(min(double(crop) * (255.0 / bg_level), 255));
            else
                norm_crop = crop;
            end

            accum_min = min(accum_min, norm_crop);
            counter   = counter + 1;
        end

        clear vid;
        fprintf('Processed fly %d: added %d break frames.\n', fly_id, length(break_frames));

    catch ME
        fprintf('Error on fly %d: %s\n', fly_id, ME.message);
    end
end

% --- Final save ---
if counter > 0
    output_path = fullfile(OUTPUT_DIR, 'global_break_posture_stack.png');
    final_img   = normalize_minmax(accum_min);
    imwrite(final_img, output_path);
    fprintf('\nSuccess! Stacked %d frames into %s\n', counter, output_path);
else
    fprintf('\nFailed: no frames were processed successfully.\n');
end
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
