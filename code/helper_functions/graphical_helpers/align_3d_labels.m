%% Script Name:
% align_3d_labels
%
% align_3d_labels(ax, pad)
%
% Rotates the x and y labels of a 3D axes so they lie along their axis edges,
% and parks them just outside the box. Call it after limits, ticks and view
% are set. Does nothing useful on a flat view such as view(90, 90).
%
%   ax    axes handle
%   pad   distance of the label from the axis, as a fraction of the other
%         axis range                                           (default: 0.18)

function align_3d_labels(ax, pad)

if nargin < 2
    pad = 0.18;
end

[az, el] = view(ax);
xl = ax.XLim; yl = ax.YLim; zl = ax.ZLim;
pb = ax.PlotBoxAspectRatio;

% The camera side decides which edges of the box are in front
cam = [sin(deg2rad(az)), -cos(deg2rad(az))];
y_edge = yl(1 + (cam(2) > 0));   % front edge in y, where the x axis sits
x_edge = xl(1 + (cam(1) > 0));   % front edge in x, where the y axis sits
z0 = zl(1);

V = viewmtx(az, el);
lo = [xl(1) yl(1) zl(1)];
rng = [diff(xl) diff(yl) diff(zl)];

xa = edge_angle(V, lo, rng, pb, [xl(1) y_edge z0], [xl(2) y_edge z0]);
ya = edge_angle(V, lo, rng, pb, [x_edge yl(1) z0], [x_edge yl(2) z0]);

dy = sign(y_edge - mean(yl)) * pad * diff(yl);
dx = sign(x_edge - mean(xl)) * pad * diff(xl);

set(ax.XLabel, 'Rotation', xa, 'HorizontalAlignment', 'center', ...
    'Position', [mean(xl), y_edge + dy, z0]);
set(ax.YLabel, 'Rotation', ya, 'HorizontalAlignment', 'center', ...
    'Position', [x_edge + dx, mean(yl), z0]);

end

function a = edge_angle(V, lo, rng, pb, p1, p2)
% Screen angle of the edge p1 -> p2, in the normalised box MATLAB projects

q1 = V * [((p1 - lo) ./ rng) .* pb, 1]';
q2 = V * [((p2 - lo) ./ rng) .* pb, 1]';

a = atan2d(q2(2) - q1(2), q2(1) - q1(1));
a = mod(a + 90, 180) - 90;   % keeps the text upright

end
