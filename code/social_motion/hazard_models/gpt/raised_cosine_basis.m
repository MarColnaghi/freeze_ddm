function B = raised_cosine_basis(x, nb, width_mult, anchor)
% RAISED_COSINE_BASIS  Linearly-spaced raised-cosine bumps tiling range(x).
%   width_mult : bump half-width in units of centre spacing (>1 ⇒ more overlap).
%   anchor     : optional; shift the (uniform) grid by <½ spacing so one centre
%                lands exactly on `anchor` (e.g. 0). Used for the baseline basis here
%                (anchor unused) but kept signature-compatible with the kernel code.
    if nargin < 3 || isempty(width_mult), width_mult = 1; end
    if nargin < 4, anchor = []; end
    x  = x(:);
    lo = min(x); hi = max(x);
    c  = linspace(lo, hi, nb);
    sp = (hi - lo) / (nb - 1);
    if ~isempty(anchor)
        [~, jn] = min(abs(c - anchor));
        c = c + (anchor - c(jn));
    end
    w = width_mult * sp;
    B = zeros(numel(x), nb);
    for j = 1:nb
        d = (x - c(j)) / w;
        in = abs(d) <= 1;
        B(in, j) = 0.5 * (1 + cos(pi * d(in)));
    end
end
