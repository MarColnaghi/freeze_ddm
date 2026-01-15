function [ll, all_pdfs] = compute_likelihood_tv_ddm(model_results, model, soc_mot_vector)

mu_1 = model_results.estimates_mean.mu1_sm;
pr = param_res(out.theta1(idx_bout));
out_s = ddm_pdf_from_trace([0, out.theta1(idx_bout)], soc_mot_vector, pr);

mu_2 = model_results.estimates_mean.mu2_sm;
pr = param_res(out.theta2(idx_bout));
out_l =  ddm_pdf_from_trace([0, out.theta2(idx_bout)], soc_mot_vector, pr);

pdf_ddm = out.pmix(idx_bout) * out_s.pdf + (1 - out.pmix(idx_bout)) * out_l.pdf;

pdf_ddm_survival = out.pmix(idx_bout) * out_s.survival + (1 - out.pmix(idx_bout)) * out_l.survival;

all_pdfs{idx_bout} = [pdf_ddm; pdf_ddm_survival];

if is_censored
    ll(idx_bout) = log(pdf_ddm_survival);
else
    ll(idx_bout) = log(pdf_ddm(freeze.durations));
end