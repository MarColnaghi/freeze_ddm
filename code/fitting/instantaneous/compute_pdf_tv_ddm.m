function [pdf] = compute_pdf_tv_ddm(out, points)
    % Set default NDT if missing
    if ~isfield(out, 'tndt')
        out.tndt = zeros(height(out), 1);
    end

    dt = 1/300;

    % 1. Get the raw (shifted but NOT yet truncated) outputs for both components
    pr1 = param_res(out.theta1, 'points', points, 'dt', dt);
    out_s = ddm_pdf_from_trace([0, out.theta1, out.tndt], out.mu1, pr1);

    pr2 = param_res(out.theta2, 'points', points, 'dt', dt);
    out_l = ddm_pdf_from_trace([0, out.theta2, out.tndt], out.mu2, pr2);

    % 2. Mix the SHIFTED (unnormalized) PDFs
    % We mix the raw populations before any truncation logic is applied.
    mixed_pdf_unnorm = out.pmix * out_s.pdf_shifted + (1 - out.pmix) * out_l.pdf_shifted;

    % 3. Mix the survival masses (unnormalized) at T_max
    % This is the mass that "never reached a bound" by the end of the simulation.
    mixed_surv_at_tmax_unnorm = out.pmix * out_s.survival_raw + (1 - out.pmix) * out_l.survival_raw; 

    % 4. Calculate the combined Normalization Factor (Total survival past T_trunc)
    idx_trunc = out_s.idx_trunc; 
    % Mass between T_trunc and T_max
    mass_in_window = sum(mixed_pdf_unnorm(idx_trunc:end)) * dt;
    % Total mass surviving T_trunc
    combined_norm_factor = max(1e-12, mass_in_window + mixed_surv_at_tmax_unnorm);

    % 5. Package Final Output
    pdf.t = out_s.t;
    
    % Normalized PDF: Probability of RT | RT > T_trunc
    pdf.ddm = mixed_pdf_unnorm / combined_norm_factor;
    
    % Normalized Survival: Probability of RT > T_max | RT > T_trunc
    % This is crucial for handling "omission" or "timeout" trials in your likelihood.
    pdf.survival = mixed_surv_at_tmax_unnorm / combined_norm_factor;
    
    % Store the factor just in case you need to debug the shift in proportions
    pdf.norm_factor = combined_norm_factor;
end