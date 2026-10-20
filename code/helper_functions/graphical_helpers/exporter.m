function exporter(fh, paths, filename, varargin)

opt = inputParser;
addParameter(opt, 'flag', true);
addParameter(opt, 'padding', 'figure');

parse(opt, varargin{:});

if opt.Results.flag

    fl = fullfile(paths.fig, filename);
    if strcmp(opt.Results.padding, 'figure') && endsWith(fl, '.pdf')
        % full figure canvas (exportgraphics crops to content in R2022b), so
        % figures with the same Position export at the same page size
        set(fh, 'PaperUnits', 'points', 'PaperPositionMode', 'auto');
        pos = get(fh, 'Position');
        set(fh, 'PaperSize', pos(3:4));
        col = fh.Color;
        fh.Color = 'none';              % transparent background
        set(fh, 'InvertHardcopy', 'off');
        print(fh, fl, '-dpdf', '-painters');
        fh.Color = col;
    else
        exportgraphics(fh, fl, 'BackgroundColor', 'none', 'ContentType', 'vector');
    end
end

