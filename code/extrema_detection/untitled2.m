load('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/extrema_detection/systematic_analysis/run01/mu5_theta6_snr60/sims_ac/d120_2bexp60_expkern15/struct_ac.mat')

for trials = 30:50
trial_1 = s(trials).boutlist;
data = table();
data.diffRT =  trial_1.rt_post_template - trial_1.rt_post_template(1);
data.diff_AcEv_before = trial_1.summed_motion_b4 - trial_1.summed_motion_b4(1);
data.diff_AcEv_during = trial_1.summed_motion_during - trial_1.summed_motion_during(1);
data.condition = ones(height(data), 1);
analyze_RT_with_accumulated_evidence(data(1:150,:))
end