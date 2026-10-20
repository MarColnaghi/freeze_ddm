function apply_font(fig, font_name)
% Sets the font of every text-bearing object in a figure (axes, labels, text,
% annotations, legends, colorbars, tiledlayout labels). Call it just before saving.

if nargin < 1 || isempty(fig), fig = gcf; end
if nargin < 2, font_name = 'Helvetica-Light'; end

set(findall(fig, '-property', 'FontName'), 'FontName', font_name);
