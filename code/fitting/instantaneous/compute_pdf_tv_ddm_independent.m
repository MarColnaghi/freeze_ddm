function [pdf] = compute_pdf_tv_ddm_independent(out, points)

    % This version uses pre-normalized PDFs (Truncation applied BEFORE mixing)
    if ~isfield(out, 'tndt')
        out.tndt = zeros(height(out), 1);
    end
    dt = 1/300;

    % 1. Get the outputs (Internally, ddm_pdf_from_trace applies its own normalization)
    pr1 = param_res(out.theta1, 'points', points, 'dt', dt);
    out_s = ddm_pdf_from_trace([0, out.theta1, out.tndt], out.mu1, pr1);
    
    pr2 = param_res(out.theta2, 'points', points, 'dt', dt);
    out_l = ddm_pdf_from_trace([0, out.theta2, out.tndt], out.mu2, pr2);

    % 2. Mix the ALREADY NORMALIZED components
    % Each out_x.pdf already integrates to 1.0 (over the window T_trunc to infinity)
    pdf.ddm = out.pmix * out_s.pdf + (1 - out.pmix) * out_l.pdf;
    
    % 3. Mix the survival probabilities
    % Each out_x.survival is already relative to its own T_trunc survival
    pdf.survival = out.pmix * out_s.survival + (1 - out.pmix) * out_l.survival;

    % 4. Package time vector
    pdf.t = out_s.t;
end