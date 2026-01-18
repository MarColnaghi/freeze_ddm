function ll = compute_likelihood(pdf, is_censored, freeze_duration)

if is_censored
    ll = log(pdf.survival);
else
    likelihood = interp1(pdf.t, pdf.ddm, freeze_duration);
    ll = log(likelihood);

end