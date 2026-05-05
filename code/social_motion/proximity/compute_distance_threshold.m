function results = compute_distance_threshold(P, D, params)
% P: pixel change (Nx1)
% D: min distance (Nx1)
% params: struct with optional fields: edges, nBins, epsilonFrac, etc.

%% --- Defaults & Input Handling ---
if ~isfield(params, 'epsilonFrac'), params.epsilonFrac = 0.15; end
if ~isfield(params, 'smoothSpan'), params.smoothSpan = 0.1; end
if ~isfield(params, 'farFrac'), params.farFrac = 0.2; end

% Priority: 1. params.edges, 2. params.nBins, 3. Default (50)
if isfield(params, 'edges')
    edges = params.edges;
    params.nBins = length(edges) - 1;
elseif isfield(params, 'nBins')
    edges = linspace(min(D, [], 'omitnan'), max(D, [], 'omitnan'), params.nBins + 1);
else
    params.nBins = 50;
    edges = linspace(min(D, [], 'omitnan'), max(D, [], 'omitnan'), params.nBins + 1);
end

% Remove NaNs
valid = ~isnan(P) & ~isnan(D);
P = P(valid);
D = D(valid);

%% 1. BINNING
binCenters = (edges(1:end-1) + edges(2:end))/2;
P_bin = nan(params.nBins,1);
counts = zeros(params.nBins,1);

for i = 1:params.nBins
    idx = D >= edges(i) & D < edges(i+1);
    if any(idx)
        P_bin(i) = median(P(idx)); % robust to outliers
        counts(i) = sum(idx);
    end
end

%% 2. SMOOTHING (LOESS)
validBins = ~isnan(P_bin);
D_fit = binCenters(validBins);
P_fit = P_bin(validBins);
P_smooth = smooth(D_fit, P_fit, params.smoothSpan, 'loess');

%% 3. BASELINE (far distances)
[~, sortIdx] = sort(D);
nFar = round(params.farFrac * length(D));
farIdx = sortIdx(end-nFar+1:end);
P_far = median(P(farIdx));

%% 4. THRESHOLD DETECTION (Minimum-based)

% find minimum (clean region)
[P_min, idx_min] = min(P_smooth);

% tolerance around minimum
epsilon = params.epsilonFrac * abs(P_min);

% find first point (from left) entering clean region
idx_thresh = find(P_smooth <= (P_min + epsilon), 1, 'first');

if isempty(idx_thresh)
    D_thresh = NaN;
else
    D_thresh = D_fit(idx_thresh);
end



%% 5. OUTPUT
results.edges = edges;
results.binCenters = binCenters;
results.P_bin = P_bin;
results.D_fit = D_fit;
results.P_smooth = P_smooth;
results.P_far = P_far;
results.threshold = D_thresh;
results.epsilon = epsilon;
results.counts = counts;
end