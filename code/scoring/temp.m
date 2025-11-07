% Assumptions:
%   - Table is T
%   - Fly ID column: T.fly
%   - Loom identifier per window: T.nloom        % (use your actual column; e.g., nloom or nloom_loomwin)
%   - Bout timing: T.onsets                      % (for ordering within a loom)
%   - Bout duration: T.durations_s               % (use T.durations if you prefer frames)

T = bouts_proc;

% 1) Sort so bouts are ordered within each loom window
T = sortrows(T, {'fly','nloom','onsets'});

% 2) Group by fly & loom
[G, Keys] = findgroups(T(:, {'fly','nloom'}));

% 3) Gather durations per (fly, loom) as a row vector (ordered by onsets thanks to the sort)
dur_cells = splitapply(@(d){d(:).'}, T.durations_s, G);

% 4) Count how many bouts happened in each (fly, loom)
bout_counts = splitapply(@numel, T.durations_s, G);

% 5A) Keep all multi-bout looms (>=2 bouts)
idx_multi = bout_counts >= 2;
PairsAll           = Keys(idx_multi, :);
PairsAll.n_bouts   = bout_counts(idx_multi);
PairsAll.durations = dur_cells(idx_multi);   % cell array: each cell is [d1 d2 ...] for that loom

figure
hold on
scatter(1, PairsAll.durations)
scatter(2, PairsAll.durations(:,2))

