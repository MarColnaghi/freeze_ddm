% Audit: do the surrounding LC6 (immobile) flies jump at the first looms?
% Expected number of jumping neighbours per recording = LC6 count in the name (4 - NorpA).
% Read-only: reports discrepancies, does not modify any data.

raw_dir = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_using_raw_files';
out_dir = '/Users/marcocolnaghi/PhD/freeze_ddm/datasets/CS';

n_looms   = 7;     % looms checked
win       = 120;   % frames after loom onset in which a jump is counted
vel_jump  = 75;    % mm/s, velocity criterion (README rule for the jumps flag)

f1 = dir(fullfile(raw_dir, '*.fly1_scored_frames.csv'));
recs = erase({f1.name}, '.fly1_scored_frames.csv')';
nrec = numel(recs);

norpa = nan(nrec, 1); lc6 = nan(nrec, 1);
for r = 1:nrec
    tok = regexp(recs{r}, '(\d)NorpA(\d)LC6ChR', 'tokens', 'once');
    norpa(r) = str2double(tok{1}); lc6(r) = str2double(tok{2});
end

jump_ev  = false(nrec, 4, n_looms);   % neighbour has jumps>0 in window after loom l
vel_ev   = false(nrec, 4, n_looms);   % neighbour has velocity>=vel_jump in window
n_starts = nan(nrec, 1);
total_flags = nan(nrec, 4);           % all jump flags over the whole recording

parfor r = 1:nrec
    opts = detectImportOptions(fullfile(raw_dir, [recs{r} '.fly1_scored_frames.csv']));
    opts = setvartype(opts, 'double');
    opts.SelectedVariableNames = {'looming_bout'};
    loom = readmatrix(fullfile(raw_dir, [recs{r} '.fly1_scored_frames.csv']), opts);
    starts = find(diff(loom) == 1) + 1;
    n_starts(r) = numel(starts);

    je = false(4, n_looms); ve = false(4, n_looms); tf = nan(1, 4);
    for k = 1:4
        fn = fullfile(raw_dir, sprintf('%s.fly%d_scored_frames.csv', recs{r}, k + 1));
        o = detectImportOptions(fn);
        o.SelectedVariableNames = {'jumps', 'velocity'};
        o = setvartype(o, 'double');
        T = readtable(fn, o);
        j = T.jumps; v = T.velocity;
        j(isnan(j)) = 0; v(isnan(v)) = 0;
        tf(k) = sum(j > 0);
        for l = 1:min(n_looms, numel(starts))
            idx = starts(l):min(starts(l) + win, numel(j));
            je(k, l) = any(j(idx) > 0);
            ve(k, l) = any(v(idx) >= vel_jump);
        end
    end
    jump_ev(r, :, :) = reshape(je, 1, 4, n_looms); vel_ev(r, :, :) = reshape(ve, 1, 4, n_looms); total_flags(r, :) = tf;
end

exp_n       = lc6;
n_L1        = sum(jump_ev(:, :, 1), 2);
n_any       = sum(any(jump_ev, 3), 2);
n_any_vel   = sum(any(vel_ev, 3), 2);
n_per_loom  = squeeze(sum(jump_ev, 2));          % nrec x n_looms
looms_per_nb = squeeze(sum(jump_ev, 3));         % nrec x 4: looms (of 7) each neighbour jumped at
match       = n_any == exp_n;
match_L1    = n_L1 == exp_n;

audit = table(recs, norpa, lc6, exp_n, n_L1, n_any, n_any_vel, match, match_L1, n_starts, ...
    n_per_loom, looms_per_nb, total_flags, ...
    'VariableNames', {'rec','norpa','lc6','expected','n_jump_loom1','n_jump_any1to7','n_vel_any1to7', ...
    'match_any','match_loom1','n_looms_found','n_per_loom','looms_per_neighbour','flags_total'});

save(fullfile(out_dir, 'lc6_jumper_audit.mat'), 'audit', 'jump_ev', 'vel_ev');
writetable(audit(:, {'rec','norpa','lc6','expected','n_jump_loom1','n_jump_any1to7','n_vel_any1to7','match_any','match_loom1','n_looms_found'}), ...
    fullfile(out_dir, 'lc6_jumper_audit.csv'));

fprintf('\nRecordings: %d, looms found != 20 in %d\n', nrec, sum(n_starts ~= 20));
for g = 0:4
    m = norpa == g;
    fprintf('\nNorpA=%d (expected %d jumpers), n=%d | exact match (any loom 1-7): %d (%.1f%%), loom 1 only: %d (%.1f%%)\n', ...
        g, 4 - g, sum(m), sum(match(m)), 100 * mean(match(m)), sum(match_L1(m)), 100 * mean(match_L1(m)));
    fprintf('  observed jumpers (any 1-7), count of recordings: ');
    for c = 0:4, fprintf('%d:%d  ', c, sum(m & n_any == c)); end
    fprintf('\n  observed jumpers (loom 1)  : ');
    for c = 0:4, fprintf('%d:%d  ', c, sum(m & n_L1 == c)); end
    fprintf('\n  velocity-criterion (any 1-7): ');
    for c = 0:4, fprintf('%d:%d  ', c, sum(m & n_any_vel == c)); end
    fprintf('\n  mean #jumper flies per loom 1..7: %s\n', mat2str(mean(n_per_loom(m, :), 1), 3));
end
