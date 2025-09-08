function apply_generic(ax, fontsize)

if nargin < 2
    fontsize = 24;
end
ax.XAxis.FontSize = fontsize;
ax.ZAxis.FontSize = fontsize;
ax.YAxis.FontSize = fontsize;
set(ax,'linewidth', 3, 'TickDir','out');
set(ax,'box','off')
set(ax,'TickLength',[0.015, 0.015])