function apply_generic(ax, varargin)

opt = inputParser;
addParameter(opt, 'font_size', 24);
addParameter(opt, 'tick_length', 0.015);
addParameter(opt, 'line_width', 3);
addParameter(opt, 'no_xticks', false); 
addParameter(opt, 'no_yticks', false); 
addParameter(opt, 'no_xlabels', false); 
addParameter(opt, 'no_ylabels', false); 

addParameter(opt, 'no_x', false); 
addParameter(opt, 'no_y', false); 
addParameter(opt, 'xlims', []); 
addParameter(opt, 'ylims', []); 
addParameter(opt, 'xticks', []); 
addParameter(opt, 'yticks', []); 
addParameter(opt, 'box', 'off'); 
addParameter(opt, 'ypos', []); 
addParameter(opt, 'xpos', []); 

parse(opt, varargin{:});

fontsize = opt.Results.font_size;
tick_length = opt.Results.tick_length;
line_width = opt.Results.line_width;
no_xticks = opt.Results.no_xticks;
no_yticks = opt.Results.no_yticks;
no_xlabels = opt.Results.no_xlabels;
no_ylabels = opt.Results.no_ylabels;
no_xaxis = opt.Results.no_x;
no_yaxis = opt.Results.no_y;
xlims = opt.Results.xlims;
ylims = opt.Results.ylims;
xthicks = opt.Results.xticks;
ythicks = opt.Results.yticks;
box = opt.Results.box;
ypos = opt.Results.ypos;
xpos = opt.Results.xpos;

% ax.YAxis is a 1x2 ruler array when yyaxis is active; set() handles both
ax.XAxis.FontSize = fontsize;
ax.ZAxis.FontSize = fontsize;
set(ax.YAxis, 'FontSize', fontsize);
set(ax,'linewidth', line_width, 'TickDir','out');
set(ax,'box', box)
set(ax,'TickLength',[tick_length, tick_length])
set(ax,'Color', 'none')

if ~isempty(xlims)
    xlim(ax, xlims)
end

% with yyaxis active, ylims/yticks target the currently active side only
if ~isempty(ylims)
    ylim(ax, ylims)
end

if ~isempty(xthicks)
    xticks(ax, xthicks)
end

if ~isempty(ythicks)
    yticks(ax, ythicks)
end

if no_xticks
    xticks(ax, []);
end

if no_yticks
    set(ax.YAxis, 'TickValues', []);
end

if no_xlabels
    xticklabels(ax, {});
end

if no_ylabels
    set(ax.YAxis, 'TickLabels', {});
end

if no_xaxis
    ax.XAxis.Visible = 'off';
end

if no_yaxis
    set(ax.YAxis, 'Visible', 'off');
end

if ~isempty(ypos) && numel(ax.YAxis) == 1
    % YAxisLocation cannot be changed while yyaxis is active
    set(ax, 'YAxisLocation', ypos)
end

if ~isempty(xpos)
    set(ax, 'XAxisLocation', xpos)
end
