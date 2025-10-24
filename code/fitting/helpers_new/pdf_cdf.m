function [pdf, cdf] = pdf_cdf(models)

% Define all possible pdf and cdf functions
all_pdf.ddm = @(ts, mu, theta, ndt) theta ./ sqrt(2 * pi * (ts - ndt) .^ 3) .* ...
    exp(-(mu .* (ts - ndt) - theta) .^ 2 ./ (2 * (ts - ndt)));

all_cdf.ddm = @(ts, mu, theta, ndt) 1 - normcdf((theta - mu .* (ts - ndt)) ./ sqrt(ts - ndt)) + ...
    exp(2 .* mu .* theta) .* normcdf((-theta - mu .* (ts - ndt)) ./ sqrt(ts - ndt));
        
all_pdf.exp = @(ts, lambda) lambda .* exp(-lambda .* ts);
all_cdf.exp = @(ts, lambda) 1 - exp(-lambda .* ts);

all_pdf.kde_pdf_interp = @(z) interp1(extra.xkde, extra.fkde, z, 'linear', 0);
all_pdf.kde_cdf_interp = @(z) interp1(extra.xkde, extra.Fkde, z, 'linear', 0);

all_pdf.ed = @(ts, theta, mu, ndt, fs) ed_vectorized_trials_log_continuous(ts - ndt, theta, mu, fs, 'pdf');
all_cdf.ed = @(ts, theta, mu, ndt, fs) ed_vectorized_trials_log_continuous(ts - ndt, theta, mu, fs, 'cdf');

% Initialize empty structs
pdf = struct();
cdf = struct();

% Validate input and populate only the requested fields
if iscell(models)
    for i = 1:length(models)
        model = models{i};
        if isfield(all_pdf, model)
            pdf.(model) = all_pdf.(model);
            cdf.(model) = all_cdf.(model);
        else
            error('Model "%s" is not recognized.', model);
        end
    end
else
    error('Input must be a cell array of model names, e.g., {''ddm'', ''exp''}.');
end
end
