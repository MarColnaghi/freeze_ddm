clear all
close all
idx_seed = randi(100);

exporting = true;
col = cmapper();
sim_params.rng = idx_seed;
rng(sim_params.rng);

% Model
model = 'dddm2';
select_run = 'run03';
gen_data = 'fr';

% Load important files
paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
load(fullfile(paths.results, 'surrogate.mat'))
bouts_proc = surrogate;

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

% 
% 
% %%
% 
% 
% for idx_fly = 1
% figure
% 
%     sm = motion_cache((idx_fly));
%     [pks, ind] = findpeaks(sm, "MinPeakHeight", 2);
% 
%     plot(sm)
%     scatter(ind, pks)
%     ylim([0 10])
%     
%     bouts_fly = bouts_proc(bouts_proc.fly == idx_fly, :);
% 
%     for idx_bout = 1:height(bouts_fly)
%         sm_bout = sm(bouts_fly.onsets(idx_bout):bouts_fly.onsets(idx_bout) + bouts_fly.durations(idx_bout));
%         if any(ind >= bouts_fly.onsets(idx_bout) & ind < bouts_fly.onsets(idx_bout) + bouts_fly.durations(idx_bout))
%             figure
%             hold on
%             plot(bouts_fly.onsets(idx_bout):bouts_fly.onsets(idx_bout) + bouts_fly.durations(idx_bout), sm_bout, 'b')
%             scatter(ind(ind >= bouts_fly.onsets(idx_bout) & ind < bouts_fly.onsets(idx_bout) + bouts_fly.durations(idx_bout)), pks(ind >= bouts_fly.onsets(idx_bout) & ind < bouts_fly.onsets(idx_bout) + bouts_fly.durations(idx_bout)))
% 
%         end
%     end
% 
% end
% 
% %%

bouts_proc = bouts_proc(bouts_proc.durations > 60,:);
d = 60;

for idx_bout = 1 % height(bouts_proc)

    sm = motion_cache(bouts_proc.fly(idx_bout));
    sm_bout = sm(bouts_proc.onsets(idx_bout):bouts_proc.onsets(idx_bout) + bouts_proc.durations(idx_bout));
    offset = length(sm_bout); onset = offset - d;
    v1 = sm_bout(onset:offset);

    sm_cell = extract_sm_cellarray(bouts_proc, motion_cache);

    for idx_comparison = 1:height(bouts_proc)
        sm_idx = sm_cell{idx_comparison};
        v2 = sm_idx(onset:offset);

        similarity(idx_comparison) = norm(v1 - v2);

    end

    figure
    hold on
    histogram(similarity, 0:0.25:50)

end




%%
figure 
hold on
for n_mov = 0:4 
    bouts_proc = surrogate(surrogate.moving_flies == n_mov, :);
    i = 0;
    for idx_bout = 1:height(bouts_proc)
        i = i + 1;
        sm = motion_cache(bouts_proc.fly(idx_bout));
        sm_bout(i, :) = sm(bouts_proc.onsets(idx_bout):bouts_proc.onsets(idx_bout) + bouts_proc.durations(idx_bout));
    end
    plot(mean(sm_bout))

end

