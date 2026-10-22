function exporter(fh, paths, filename, varargin)

opt = inputParser;
addParameter(opt, 'flag', true);
addParameter(opt, 'padding', 'figure');
addParameter(opt, 'font_name', 'Helvetica-Light');

parse(opt, varargin{:});

if opt.Results.flag

    apply_font(fh, opt.Results.font_name);

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
    
    elseif strcmp(opt.Results.padding, 'figure') && endsWith(fl, '.svg')
        % full figure canvas (exportgraphics crops to content in R2022b), so
        % figures with the same Position export at the same page size
        set(fh, 'PaperUnits', 'points', 'PaperPositionMode', 'auto');
        pos = get(fh, 'Position');
        set(fh, 'PaperSize', pos(3:4));
        col = fh.Color;
        fh.Color = 'none';              % transparent background
        set(fh, 'InvertHardcopy', 'off');
        print(fh, fl, '-dsvg', '-painters');
        fh.Color = col;
        % print writes font-family:'Helvetica' whatever FontName is, losing the
        % weight; put the requested face back so Illustrator picks it up
        txt = fileread(fl);
        txt = strrep(txt, "font-family:'Helvetica'", ...
            sprintf("font-family:'%s'", strrep(opt.Results.font_name, '-', ' ')));
        fid = fopen(fl, 'w');
        fwrite(fid, txt);
        fclose(fid);
    else
        exportgraphics(fh, fl, 'BackgroundColor', 'none', 'ContentType', 'vector');
    end
end

