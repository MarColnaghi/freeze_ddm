function results = zero_pixel_distance_analysis(P, D)

% remove NaNs
valid = ~isnan(P) & ~isnan(D);
P = P(valid);
D = D(valid);

% find exact zeros
zero_idx = (P == 0);

D_zero = D(zero_idx);

% basic stats
results.count = sum(zero_idx);
results.fraction = mean(zero_idx);

if isempty(D_zero)
    warning('No zero pixel change frames found.');
    results.D_zero = [];
    results.D_min = NaN;
    results.D_p5 = NaN;
    return;
end

% distances when pixel change is zero
results.D_zero = D_zero;

% robust lower bound
results.D_min = min(D_zero);
results.D_p5 = prctile(D_zero, 5);   % better than min