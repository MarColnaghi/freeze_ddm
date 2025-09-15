function sm_cell = extract_sm_cellarray(bouts, motion_cache)

for idx_bout = 1:height(bouts)
    sm_vector = motion_cache(bouts.fly(idx_bout));
    ons = bouts.onsets(idx_bout);
    off = bouts.onsets(idx_bout) + round(bouts.durations_s(idx_bout) * 60);
    sm_cell{idx_bout} = sm_vector(ons:off);
end
