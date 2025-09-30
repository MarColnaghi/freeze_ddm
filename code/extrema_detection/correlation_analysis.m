%calculate_corrs()

n = 60;

%
for idx_freezes = 1:size(s,2)
    temp_s = s(idx_freezes).boutlist;
    similarities(idx_freezes) = temp_s.closest_similarity(n);
    corrmat = corrcoef(temp_s.summed_motion_b4(1:n), temp_s.rt_post_template(1:n));
    correlations(idx_freezes) = corrmat(2);
end

figure

histogram(correlations, -0.75:0.05:0.75)

%%
for idx_freezes = 1:15
    temp_s = s(idx_freezes).boutlist;
    similarities(idx_freezes) = temp_s.closest_similarity(n);
    corrmat = corrcoef(temp_s.summed_motion_b4(1:n), temp_s.rt_post_template(1:n));
    correlations(idx_freezes) = corrmat(2);
    figure
    scatter(temp_s.summed_motion_b4(1:n), temp_s.rt_post_template(1:n))
end

