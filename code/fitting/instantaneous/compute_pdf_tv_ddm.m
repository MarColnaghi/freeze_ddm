function [pdf] = compute_pdf_tv_ddm(out)

pr = param_res(out.theta1);
out_s = ddm_pdf_from_trace([0, out.theta1 out.tndt], out.mu1, pr);

pr = param_res(out.theta2);
out_l =  ddm_pdf_from_trace([0, out.theta2 out.tndt], out.mu2, pr);

pdf.ddm = out.pmix * out_s.pdf + (1 - out.pmix) * out_l.pdf;

pdf.survival = out.pmix * out_s.survival + (1 - out.pmix) * out_l.survival;

pdf.t = out_s.t;