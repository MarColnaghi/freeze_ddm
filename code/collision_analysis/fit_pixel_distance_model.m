function model = fit_pixel_distance_model(P, D)

% remove NaNs
valid = ~isnan(P) & ~isnan(D);
P = P(valid);
D = D(valid);

% model: a*exp(-D/lambda) + b + c*D
model_fun = @(params, D) ...
    params(1) * exp(-D / params(2)) + params(3) + params(4)*D;

% initial guesses (important!)
a0 = max(P) - min(P);   % amplitude
lambda0 = median(D);    % decay scale
b0 = min(P);            % baseline
c0 = 0;                 % slope

params0 = [a0, lambda0, b0, c0];

% fit
opts = optimset('Display','off');

params_fit = fminsearch(@(p) ...
    sum((P - model_fun(p, D)).^2), params0, opts);

% store results
model.params = params_fit;
model.fun = model_fun;

end