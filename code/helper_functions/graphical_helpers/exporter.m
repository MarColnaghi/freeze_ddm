function exporter(fh, paths, filename, varargin)

opt = inputParser;
addParameter(opt, 'flag', true);

parse(opt, varargin{:});

if opt.Results.flag

    fl = fullfile(paths.fig, filename);
    exportgraphics(fh, fl, 'BackgroundColor', 'none', 'ContentType', 'vector');
end

