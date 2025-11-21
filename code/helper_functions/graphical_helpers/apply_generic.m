function apply_generic(ax, varargin)

opt = inputParser;
addParameter(opt, 'font_size', 24);
addParameter(opt, 'tick_length', 0.015);
addParameter(opt, 'line_width', 3);
addParameter(opt, 'no_xticks', false); 
addParameter(opt, 'no_yticks', false); 
addParameter(opt, 'no_x', false); 
addParameter(opt, 'no_y', false); 
addParameter(opt, 'xlims', []); 
addParameter(opt, 'ylims', []); 

parse(opt, varargin{:});

fontsize = opt.Results.font_size;
tick_length = opt.Results.tick_length;
line_width = opt.Results.line_width;
no_xticks = opt.Results.no_xticks;
no_yticks = opt.Results.no_yticks;
no_xaxis = opt.Results.no_x;
no_yaxis = opt.Results.no_y;
xlims = opt.Results.xlims;
ylims = opt.Results.ylims;

ax.XAxis.FontSize = fontsize;
ax.ZAxis.FontSize = fontsize;
ax.YAxis.FontSize = fontsize;
set(ax,'linewidth', line_width, 'TickDir','out');
set(ax,'box','off')
set(ax,'TickLength',[tick_length, tick_length])

if ~isempty(xlims)
    xlim(xlims)
end

if ~isempty(ylims)
    ylim(ylims)
end

if no_xticks
    xticks([]);
end

if no_yticks
    yticks([]);
end

if no_xaxis
    ax.XAxis.Visible = 'off';
end

if no_yaxis
    ax.YAxis.Visible = 'off';
end