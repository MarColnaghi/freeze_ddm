function mat = cache2mat(cache, varargin)

opt = inputParser;
addParameter(opt, 'selected_flies', 1:size(cache, 1));
parse(opt, varargin{:});

selected_flies = opt.Results.selected_flies;

vals = values(cache, num2cell(selected_flies));
mat  = horzcat(vals{:})';