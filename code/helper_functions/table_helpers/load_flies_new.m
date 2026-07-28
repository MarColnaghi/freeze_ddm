function [bouts, soc_mot] = load_flies_new(thresholds, varargin)


opt = inputParser;
addParameter(opt, 'save', 'both');
addParameter(opt, 'paths', '');
addParameter(opt, 'imfirst', false);

parse(opt, varargin{:});

%%

% Each file is named based on the genotype of the focal fly (1CS) and of the surrounding flies (number of moving NorpA flies and number of freezing flies LC6Chr flies),
% as well as the loom speed they were exposed to (25cm or 50cm). Within each file you will find the following columns: 1. Frame number, 2. Looming bout (1 every time
% the loom is on - should have 20 x 30 frames, which is the loom duration), 3. Velocity (actually the speed of the focal fly), 4. Walk_bout (0 or 1), 5. Freeze_bout (0 or 1),
% 5. Low_vel (low velocity behaviours 0 or 1), 6. Sum_motion (the summed motion cue the other surrounding, controlled, flies produce.

%% Output columns
%
% One row per behavioural bout, per fly. Each file's bouts are built up as `z`
% in the loop below and concatenated, so the order here is assignment order.
%
%   fly                 fly ID, parsed from the .flyN.csv suffix
%   type                1 immobility, 0 mobility
%   le                  placeholder, filled in by bouts_formatting
%   period              0 baseline, 1 loom block
%   onsets_loomaligned  bout onset relative to the previous loom
%   durations           bout duration in frames
%   nloom               index of the previous loom, 1-20 within its block
%   sloom               loom speed, 25 or 50 cm
%   bout_with_loom      number of looms starting inside the bout
%   frozen_start        per-loom flag, broadcast to every bout sharing that
%                       nloom: true when a loom started an immobility bout of
%                       at least 30 frames
%   avg_sm              mean sum_motion over the bout, capped at 630 frames
%   avg_ss              mean sum_speed over the bout, capped at 630 frames
%   avg_sa              mean sum_angle over the bout, capped at 630 frames
%   avg_fs_1s           mean focal velocity over the second before bout onset,
%                       NaN when the bout starts before frame 61
%   avg_fs_loom         mean focal velocity over the second before the previous
%                       loom, NaN when that loom lands before frame 61
%   avg_fs              mean focal velocity over the whole bout, uncapped
%   avg_pc              mean pixelchange over the whole bout, uncapped
%   genotype            1 CS. 2 LC11TNT and 3 EmptyTNT come from the dataset 2
%                       branch, which is unreachable under `for dataset = 0`
%   moving_flies        number of moving flies around the focal fly, 0-4
%   onsets              bout onset in frames
%   ends                one past the bout's last frame, so a bout spans
%                       onsets : ends - 1
%   onsets_loomwin      bout onset relative to the next loom, rewritten by
%                       bouts_formatting
%   nloom_loomwin       loom index aligned to the next loom, likewise
%   loom_ts             frame of the previous loom
%   loom_ts_n           frame of the next loom
%   loom_durs           span of the previous loom counted in frame gaps rather
%                       than frames, so one short. Nothing reads it
%   jump_flies          number of surrounding flies jumping at bout onset
%
% bouts_formatting then adds `id` and overwrites `le`, `onsets_loomwin` and
% `nloom_loomwin`.
%
% The 630-frame cap on avg_sm / avg_ss / avg_sa is deliberate, not an
% oversight. Those three are the social variables, describing what the
% surrounding flies did, and 630 frames is about one loom cycle at 60 fps: a
% ~570 frame ITI plus the 30 frame loom. It is the same window as `fine` in
% impose_contact_threshold. avg_fs and avg_pc describe the focal fly's own
% behaviour during the bout and so run the full length. The consequence worth
% remembering is that for any bout longer than 630 frames the two groups are
% averaged over different supports, so they are not directly comparable.

soc_mot = table();
soc_mot.ts_sm = zeros(0);

% Restored on the way out, including on error, so a failed run does not leave
% the warning suppressed for the rest of the session
prev_warn = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames'); % readtable warning
restore_warn = onCleanup(@() warning(prev_warn)); %#ok<NASGU>

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

for dataset = 0

    % Every dir and readtable call below prefixes `directory`, so there is no
    % need to cd in here and leave the session somewhere else on the way out
    directory = sprintf('/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_%d/', dataset);
    n_moving = 0:4;

    if dataset == 0 
        token = '1CS%dNorpA%dLC6ChR_5F-%dcm*';
        speeds = [25,50];
        genotype = 1;

    elseif dataset == 2
        token = '*%dTrh_%dNorpA-%dcm*';
        speeds = 50;

    end

    for idx_sloom = 1:length(speeds)

        for moving_flies = 1:length(n_moving)

            if dataset == 0
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
                loom_starts = loom_frames([true; diff(loom_frames)>1])+1;
                loom_ends = loom_frames([diff(loom_frames)>1;true])+1;

                % The two masks used to disagree, >2 against >1, so a single
                % dropped frame mid-loom surfaced as a dimension mismatch on
                % the subtraction below. Check the structure directly instead
                assert(numel(loom_starts) == 20 && numel(loom_ends) == 20, ...
                    'load_flies_new:loomCount', ...
                    'Expected 20 looms in %s, found %d starts and %d ends.', ...
                    listing(kf).name, numel(loom_starts), numel(loom_ends));

                loom_durs = repmat(loom_ends - loom_starts, 2, 1);

                event_indices = repmat(loom_starts, 2, 1);
                event_indices(1:20,:) = event_indices(1:20,:) - 18001;

                % interp1 needs a monotonic grid. The -18001 shift on the first
                % block only keeps it monotonic while the baseline is longer
                % than the span of the loom block
                assert(all(diff(event_indices) > 0), 'load_flies_new:eventOrder', ...
                    'Loom event grid is not strictly increasing in %s.', listing(kf).name);

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

                % A bout with no previous loom lands here as nloom == 0. That
                % is then dropped by the histcounts edges further down and only
                % surfaces as a height mismatch on the loom_durs assignment, a
                % long way from the cause
                assert(all(nloom > 0), 'load_flies_new:onsetBeforeGrid', ...
                    '%d bouts in %s have no preceding loom event (earliest onset %d, first event %d).', ...
                    sum(nloom == 0), listing(kf).name, min(onsets), event_indices(1));

                z.period = nloom > 20;

                z.onsets_loomaligned =  onsets - loom_ts_previous;
                z.durations = run_lengths;

                z.nloom = nloom;
                z.sloom = loom_speed .* ones(length(run_lengths), 1);

                % Bounds match the bout's actual span, onsets : run_ends - 1.
                % The old form paired a strict > against onsets with a strict <
                % against run_ends, so a loom landing on the bout's first frame
                % was dropped while one on its last frame was kept
                z.bout_with_loom = sum((event_indices >= onsets') & (event_indices <= (run_ends - 1)'), 1)';

                if any(z.bout_with_loom > 0 & z.type == 1 & z.durations >= 30)
                    frozen_starts = unique(z.nloom(z.bout_with_loom > 0 & z.type == 1 & z.durations >= 30) + z.bout_with_loom(z.bout_with_loom > 0 & z.type == 1 & z.durations >= 30));
                    % long_freezes = z.bout_with_loom == 1 & z.type == 1 & z.durations >= 30;
                    z.frozen_start = ismember(z.nloom, frozen_starts);
                else
                    z.frozen_start = false(size(z.nloom)); % Ensure correct output format
                end

                capped_lengths = min(run_lengths, 630);
                z.avg_sm = compute_means(Fly1.sum_motion(1:end), onsets , onsets + capped_lengths - 1);
                z.avg_ss = compute_means(Fly1.sum_speed(1:end), onsets, onsets + capped_lengths - 1);
                z.avg_sa = compute_means(Fly1.sum_angle(1:end), onsets, onsets + capped_lengths - 1);

                l.ts_sm = arrayfun(@(start_idx, end_idx) ...
                    Fly1.sum_motion(start_idx:end_idx), ...
                    run_ends - run_lengths, run_ends - 1, ...
                    'UniformOutput', false);

                z.avg_fs_1s = nan(length(run_lengths), 1);
                z.avg_fs_1s(onsets >= 61, :) = compute_means(Fly1.velocity(1:end), onsets(onsets >= 61, :) - 60, onsets(onsets >= 61, :) - 1);

                % Needs a full second of history before the loom, so the guard
                % is >= 61 to match avg_fs_1s above. The old > 0 admitted looms
                % at frames 1-60, where the window start goes non-positive and
                % the indexing throws
                z.avg_fs_loom = nan(length(run_lengths), 1);
                has_pre_loom = loom_ts_previous >= 61;
                z.avg_fs_loom(has_pre_loom) = compute_means(Fly1.velocity(1:end), loom_ts_previous(has_pre_loom) - 60, loom_ts_previous(has_pre_loom) - 1);

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

                jumps = nan(4, height(Fly1));

                for idx_surr_fly = 1:4
                    jumps(idx_surr_fly, :) = movmean(Fly1.(sprintf('speed_sur_fly_%d',idx_surr_fly)), 10) > 20;
                end

                z.jump_flies = interp1(1:height(Fly1), sum(jumps), onsets);

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

% Guarding on nargin was wrong: with varargin it counts every name-value
% element, so passing any option at all sent this block into a field reference
% on the default '' and threw. Test for the thing actually needed instead.
% Saving by full path also removes the second cd
if ~isempty(opt.Results.paths)
    bouts_file   = fullfile(opt.Results.paths.dataset, 'bouts.mat');
    soc_mot_file = fullfile(opt.Results.paths.dataset, 'soc_mot.mat');

    switch opt.Results.save
        case 'both'
            save(bouts_file, 'bouts', '-v7.3')
            save(soc_mot_file, 'soc_mot', '-v7.3')
        case 'soc_mot'
            save(soc_mot_file, 'soc_mot', '-v7.3')   % used to save `bouts` here
        case 'bouts'
            save(bouts_file, 'bouts', '-v7.3')
    end
end