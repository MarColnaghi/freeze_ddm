function ll = compute_likelihood(pdf, is_censored, duration)

if is_censored
    ll = log(pdf.survival);
else
    likelihood = interp1(pdf.t, pdf.ddm, duration);
    ll = log(likelihood);

end