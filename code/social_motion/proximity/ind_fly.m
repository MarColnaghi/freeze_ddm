%% ===============================
% INPUTS (you must have these)
% P       : pixel change (Nx1)
% D       : distance (Nx1)
% fly_id  : fly index (Nx1)
%% ===============================

%% PARAMETERS
nBins = 40;
farFrac = 0.2;
qTail = 95;        % percentile for defining contamination
pThresh = 0.05;    % probability threshold
nBoot = 50;
nShuffle = 50;

flies = unique(fly_id);
nFlies = length(flies);

D_thresh = nan(nFlies,1);
D_thresh_boot = nan(nFlies,nBoot);
D_thresh_null = nan(nFlies,nShuffle);

all_curves = cell(nFlies,1);

%% ===============================
% MAIN LOOP (per fly)
%% ===============================

for i = 1:nFlies

    idx = fly_id == flies(i);
    P_i = P(idx,:);
    D_i = D(idx,:);

    % Remove NaNs
    valid = ~isnan(P_i) & ~isnan(D_i);
    P_i = P_i(valid);
    D_i = D_i(valid);

    if length(P_i) < 50
        continue;
    end

    %% ---- BASELINE (far distances) ----
    [~, sortIdx] = sort(D_i);
    nFar = round(farFrac * length(D_i));
    farIdx = sortIdx(end-nFar+1:end);

    P_far = P_i(farIdx);

    % Define contamination threshold from far distances
    thr = prctile(P_far, qTail);

    %% ---- EXCEEDANCE ----
    isContam = P_i > thr;

    %% ---- BINNING ----
    edges = linspace(min(D_i), max(D_i), nBins+1);
    binCenters = (edges(1:end-1)+edges(2:end))/2;

    p_bin = nan(nBins,1);

    for b = 1:nBins
        idx_bin = D_i >= edges(b) & D_i < edges(b+1);
        if sum(idx_bin) > 5
            p_bin(b) = mean(isContam(idx_bin));
        end
    end

    %% ---- SMOOTH ----
    validBins = ~isnan(p_bin);
    D_fit = binCenters(validBins);
    p_smooth = smooth(D_fit, p_bin(validBins), 0.2, 'loess');

    all_curves{i} = [D_fit(:), p_smooth(:)];

    %% ---- THRESHOLD ----
    idx_thresh = find(p_smooth > pThresh, 1, 'first');

    if ~isempty(idx_thresh)
        D_thresh(i) = D_fit(idx_thresh);
    end

    %% ===============================
    % BOOTSTRAP
    %% ===============================
    for b = 1:nBoot
        bootIdx = randsample(length(P_i), length(P_i), true);

        P_b = P_i(bootIdx);
        D_b = D_i(bootIdx);

        % recompute threshold quickly (no smoothing for speed)
        isContam_b = P_b > thr;

        p_bin_b = nan(nBins,1);
        for bb = 1:nBins
            idx_bin = D_b >= edges(bb) & D_b < edges(bb+1);
            if sum(idx_bin) > 5
                p_bin_b(bb) = mean(isContam_b(idx_bin));
            end
        end

        validBins = ~isnan(p_bin_b);
        D_fit_b = binCenters(validBins);
        p_smooth_b = smooth(D_fit_b, p_bin_b(validBins), 0.2, 'loess');

        idx_thresh_b = find(p_smooth_b > pThresh, 1, 'first');

        if ~isempty(idx_thresh_b)
            D_thresh_boot(i,b) = D_fit_b(idx_thresh_b);
        end
    end

    %% ===============================
    % SHUFFLE NULL
    %% ===============================
    for s = 1:nShuffle
        P_shuff = P_i(randperm(length(P_i)));

        isContam_s = P_shuff > thr;

        p_bin_s = nan(nBins,1);
        for bb = 1:nBins
            idx_bin = D_i >= edges(bb) & D_i < edges(bb+1);
            if sum(idx_bin) > 5
                p_bin_s(bb) = mean(isContam_s(idx_bin));
            end
        end

        validBins = ~isnan(p_bin_s);
        D_fit_s = binCenters(validBins);
        p_smooth_s = smooth(D_fit_s, p_bin_s(validBins), 0.2, 'loess');

        idx_thresh_s = find(p_smooth_s > pThresh, 1, 'first');

        if ~isempty(idx_thresh_s)
            D_thresh_null(i,s) = D_fit_s(idx_thresh_s);
        end
    end

end

%% ===============================
% SUMMARY STATS
%% ===============================

mean_thresh = mean(D_thresh,'omitnan');
std_thresh = std(D_thresh,'omitnan');
cv_thresh = std_thresh / mean_thresh;

fprintf('Mean threshold: %.2f\n', mean_thresh);
fprintf('Std threshold: %.2f\n', std_thresh);
fprintf('CV: %.2f\n', cv_thresh);

%% ===============================
% PLOTS
%% ===============================

% Histogram
figure;
histogram(D_thresh,20);
xlabel('Threshold distance');
ylabel('Number of flies');
title('Per-fly threshold distribution');

% Boxplot
figure;
boxplot(D_thresh);
ylabel('Threshold distance');
title('Threshold variability');

% Example curves
figure; hold on;
nPlot = min(10,nFlies);

for i = 1:nPlot
    if ~isempty(all_curves{i})
        plot(all_curves{i}(:,1), all_curves{i}(:,2));
    end
end

xlabel('Distance');
ylabel('P(contamination)');
title('Per-fly contamination curves');

% Bootstrap variability
figure;
histogram(D_thresh_boot(:),30);
xlabel('Bootstrapped thresholds');
title('Bootstrap distribution');

% Null comparison
figure;
histogram(D_thresh,20); hold on;
histogram(D_thresh_null(:),20);
legend('Real','Shuffled');
xlabel('Threshold distance');
title('Real vs Null thresholds');