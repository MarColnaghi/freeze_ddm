function [bouts, soc_mot] = load_flies_new(thresholds, varargin)


opt = inputParser;
addParameter(opt, 'save', 'both');
addParameter(opt, 'paths', '');
addParameter(opt, 'imfirst', false);
addParameter(opt, 'edit_filename', false);

parse(opt, varargin{:});

%%

% Each file is named based on the genotype of the focal fly (1CS) and of the surrounding flies (number of moving NorpA flies and number of freezing flies LC6Chr flies),
% as well as the loom speed they were exposed to (25cm or 50cm). Within each file you will find the following columns: 1. Frame number, 2. Looming bout (1 every time
% the loom is on - should have 20 x 30 frames, which is the loom duration), 3. Velocity (actually the speed of the focal fly), 4. Walk_bout (0 or 1), 5. Freeze_bout (0 or 1),
% 5. Low_vel (low velocity behaviours 0 or 1), 6. Sum_motion (the summed motion cue the other surrounding, controlled, flies produce.

%% Initialize big table

bouts = table();

bouts.fly = zeros(0); % fly ID
bouts.genotype = zeros(0); % genotype of the fly
bouts.moving_flies = zeros(0); % number of moving flies around focal fly

bouts.sloom = zeros(0); % loom speed in cm/s
bouts.nloom = zeros(0); % loom number
bouts.nloom_loomwin = zeros(0); % loom number aligned to loom presentation
bouts.period = zeros(0); % baseline or loom

bouts.loom_start = zeros(0); % frame of loom start
bouts.loom_dur = zeros(0); % duration of each loom in frames

bouts.onsets = zeros(0); % bout onset in frames
bouts.onsets_loomaligned = zeros(0); % bout onset in frames aligned to previous loom presentation
bouts.onsets_loomwin = zeros(0); % bout onset in frames aligned to next loom presentation
bouts.le = zeros(0);

bouts.durations = zeros(0); % bout duration in frames
bouts.ends = zeros(0); % bout offset in frames

bouts.type = zeros(0); % bout type (1 immobility, 0 mobility)
bouts.bout_with_loom = zeros(0); % bout had a loom presentation inside
bouts.frozen_start = zeros(0); % bout had a loom presentation inside

bouts.avg_sm = zeros(0); % average focal fly velocity between loom start and freeze start
bouts.avg_fs = zeros(0); % average focal fly velocity between one second before loom until loom
bouts.avg_fs_1s = zeros(0); % average focal fly velocity between one second before loom until loom
bouts.avg_pc = zeros(0); % average focal fly velocity between one second before freeze until freeze

bouts.filename = zeros(0);

soc_mot = table();
soc_mot.ts_sm = zeros(0);

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames'); %turn off readtable warning

%% Enter the directory and extract properties from

% Function to compute means
compute_means = @(data, start_idx, end_indx) ...
    cell2mat(arrayfun(@(start, ends) mean(data(start:ends)), ...
    start_idx, end_indx, 'UniformOutput', false));

extract_vec = @(data, start_idx, end_indx) ...
    cell2mat(arrayfun(@(start, ends) data(start:ends), ...
    start_idx, end_indx, 'UniformOutput', false));

bouts = table();

fly_id = 0;

for dataset = 1

    directory = sprintf('/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_%d/', dataset);
    cd(directory)
    n_moving = 0:4;

    if dataset == 1
        token = '1CS%dNorpA%dLC6ChR_5F-%dcm*';
        speeds = [25,50];
        genotype = 1;

    elseif dataset == 2
        token = '*%dTrh_%dNorpA-%dcm*';
        speeds = 50;

    end

    for idx_sloom = 1:length(speeds)

        for moving_flies = 1:length(n_moving)

            if dataset == 1
                listing = dir([directory, sprintf(token, n_moving(moving_flies), 4-n_moving(moving_flies), speeds(idx_sloom))]);

            elseif dataset == 2
                listing = dir([directory, sprintf(token, 4-n_moving(moving_flies), n_moving(moving_flies), speeds(idx_sloom))]);

            end

            loom_speed = speeds(idx_sloom);

            for kf = 1:length(listing)

                Fly1 = readtable([directory, listing(kf).name]);
                Fly1.Properties.VariableNames{1} = 'frame';
                fly_id = fly_id + 1;
                tokens = regexp(listing(kf).name, '\.fly(\d+)\.csv', 'tokens');
                fly_number = str2double(tokens{1}{1});


                if dataset == 2
                    match = regexp(listing(kf).name, '-(\w+TNT)_', 'tokens');
                    if strcmp( match{1}{1}, 'LC11TNT')
                        genotype = 2;
                    elseif strcmp( match{1}{1}, 'EmptyTNT')
                        genotype = 3;
                    end
                end

                loom_frames = Fly1.frame(Fly1.looming_bout==1);
                loom_starts = loom_frames([true;diff(loom_frames)>2])+1;
                loom_ends = loom_frames([diff(loom_frames)>1;true])+1;
                loom_durs = repmat(loom_ends - loom_starts, 2, 1);

                event_indices = repmat(loom_starts, 2, 1);
                event_indices(1:20,:) = event_indices(1:20,:) - 18001;

                imm_frames = Fly1.pixelchange < thresholds.pc;

                % Filling in is done here - might make sense potentially to
                % switch around the filling in process. used to be small
                % imm first

                if opt.Results.imfirst
                    imm_frames = bwareaopen(imm_frames, thresholds.fill_in_imm); % Remove small immobile bouts
                    imm_frames = ~bwareaopen(~imm_frames, thresholds.fill_in_mob); % Remove small moving bouts
                else
                    imm_frames = ~bwareaopen(~imm_frames, thresholds.fill_in_mob); % Remove small moving bouts
                    imm_frames = bwareaopen(imm_frames, thresholds.fill_in_imm); % Remove small immobile bouts
                end

                diff_ts = [0; diff(imm_frames)]; % Detect changes
                run_ends = [find(diff_ts ~= 0); length(imm_frames)]; % Indices of run ends
                run_lengths = diff([1; run_ends]); % Calculate run lengths
                run_values = imm_frames(run_ends - 1); % Values of the runs (0 or 1)

                z = table();
                l = table();

                z.fly = fly_number .* ones(length(run_lengths), 1);
                z.type = run_values;
                z.le = zeros(length(run_lengths), 1);

                onsets = run_ends - run_lengths;
                loom_ts_previous = interp1(event_indices, event_indices, onsets, 'previous', 'extrap');
                loom_ts_next = interp1(event_indices, event_indices, onsets, 'next', 'extrap');
                [~, nloom] = ismember(loom_ts_previous, event_indices);
                z.period = nloom > 20;

                z.onsets_loomaligned =  onsets - loom_ts_previous;
                z.durations = run_lengths;

                z.nloom = nloom;
                z.sloom = loom_speed .* ones(length(run_lengths), 1);

                z.bout_with_loom = sum((event_indices > onsets') & (event_indices < run_ends'), 1)';

                if any(z.bout_with_loom > 0 & z.type == 1 & z.durations >= 30)
                    frozen_starts = unique(z.nloom(z.bout_with_loom > 0 & z.type == 1 & z.durations >= 30) + z.bout_with_loom(z.bout_with_loom > 0 & z.type == 1 & z.durations >= 30));
                    % long_freezes = z.bout_with_loom == 1 & z.type == 1 & z.durations >= 30;
                    z.frozen_start = ismember(z.nloom, frozen_starts);
                else
                    z.frozen_start = false(size(z.nloom)); % Ensure correct output format
                end

                capped_lengths = min(run_lengths, 630);
                z.avg_sm = compute_means(Fly1.sum_motion(1:end), onsets, onsets + capped_lengths - 1);

                l.ts_sm = arrayfun(@(start_idx, end_idx) ...
                    Fly1.sum_motion(start_idx:end_idx), ...
                    run_ends - run_lengths, run_ends - 1, ...
                    'UniformOutput', false);

                z.avg_fs_1s = nan(length(run_lengths), 1);
                z.avg_fs_1s(onsets >= 61, :) = compute_means(Fly1.velocity(1:end), onsets(onsets >= 61, :) - 60, onsets(onsets >= 61, :));
                z.avg_fs = compute_means(Fly1.velocity(1:end), run_ends - run_lengths, run_ends - 1);
                z.avg_pc = compute_means(Fly1.pixelchange(1:end), run_ends - run_lengths, run_ends - 1);

                z.genotype = genotype .* ones(length(run_lengths), 1);
                z.moving_flies = n_moving(moving_flies) .* ones(length(run_lengths), 1);

                z.onsets = onsets;
                z.ends = run_ends;

                z.onsets_loomwin =  onsets - loom_ts_next;
                z.nloom_loomwin = z.nloom;

                z.loom_ts = loom_ts_previous;
                z.loom_ts_n = loom_ts_next;
                z.loom_durs = repelem(loom_durs, histcounts(nloom, 0.5:1:40.5));

                soc_mot = [soc_mot; l];
                bouts = [bouts; z];

                fprintf('i%d-m%d-pc%d. Fly %d out of 983. \n', thresholds.fill_in_imm - 1, thresholds.fill_in_mob - 1, thresholds.pc, fly_id)
            end
        end
    end
end


bouts.nloom(bouts.nloom > 20) = bouts.nloom(bouts.nloom > 20) - 20;
bouts.nloom_loomwin(bouts.nloom_loomwin > 20) = bouts.nloom_loomwin(bouts.nloom_loomwin > 20) - 20;

[bouts, soc_mot] = bouts_formatting(bouts, thresholds, soc_mot);

if nargin > 1
    cd(opt.Results.paths.dataset)
    if strcmp(opt.Results.save, 'both')
        save('bouts.mat','bouts', '-v7.3')
        save('soc_mot.mat','soc_mot', '-v7.3')
    elseif strcmp(opt.Results.save, 'soc_mot')
        save('soc_mot.mat','bouts', '-v7.3')
    elseif strcmp(opt.Results.save, 'bouts')
        save('bouts.mat','bouts', '-v7.3')
    else
    end
end