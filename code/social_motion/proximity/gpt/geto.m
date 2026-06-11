function v = geto(s, f, d)
% GETO  Fetch field f from struct s, else default d (treats [] as missing).
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
