#include "mex.h"
#include <vector>
#include <cmath>
#include <algorithm>

// --- TDMA Solver (Same optimized version) ---
void solve_tdma(int n, const double* a, const double* b, const double* c, 
                const std::vector<double>& d, std::vector<double>& x,
                std::vector<double>& cp, std::vector<double>& dp) {
    double* cp_ptr = cp.data();
    double* dp_ptr = dp.data();
    const double* d_ptr = d.data();

    double b0_inv = 1.0 / b[0];
    cp_ptr[0] = c[0] * b0_inv;
    dp_ptr[0] = d_ptr[0] * b0_inv;

    for (int i = 1; i < n - 1; ++i) {
        double temp = 1.0 / (b[i] - a[i - 1] * cp_ptr[i - 1]);
        cp_ptr[i] = c[i] * temp;
        dp_ptr[i] = (d_ptr[i] - a[i - 1] * dp_ptr[i - 1]) * temp;
    }
    
    int last = n - 1;
    double temp = 1.0 / (b[last] - a[last - 1] * cp_ptr[last - 1]);
    dp_ptr[last] = (d_ptr[last] - a[last - 1] * dp_ptr[last - 1]) * temp;

    double* x_ptr = x.data();
    x_ptr[last] = dp_ptr[last];
    for (int i = n - 2; i >= 0; --i) {
        x_ptr[i] = dp_ptr[i] - cp_ptr[i] * x_ptr[i + 1];
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) mexErrMsgIdAndTxt("leaky_pde_robust:invalidNumInputs", "Inputs: drifts, params");

    double* drifts = mxGetPr(prhs[0]);
    size_t n_steps = mxGetNumberOfElements(prhs[0]);
    
    double* params = mxGetPr(prhs[1]);
    double dt = params[0];
    double dx = params[1];
    double sigma_w_sq = params[2];
    double x_min = params[4];
    double grid_size_dbl = params[5];
    int grid_size = (int)grid_size_dbl;
    int start_idx = (int)params[6] - 1;
    
    double lambda = 0.0;
    if (mxGetNumberOfElements(prhs[1]) > 7) lambda = params[7];

    // --- OUTPUT 1: The PDF (Likelihood over time) ---
    plhs[0] = mxCreateDoubleMatrix(n_steps, 1, mxREAL);
    double* likelihood = mxGetPr(plhs[0]);
    
    // --- OUTPUT 2: The Survival Probability (Scalar) ---
    // This is crucial for censoring
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* survival_out = mxGetPr(plhs[1]);

    double C_diff = (sigma_w_sq / 4.0) * (dt / (dx * dx));
    double drift_mult = dt / (4.0 * dx); 

    // Memory
    std::vector<double> p_storage(grid_size, 0.0);
    std::vector<double> p_new_storage(grid_size);
    std::vector<double>* p = &p_storage;
    std::vector<double>* p_new = &p_new_storage;
    (*p)[start_idx] = 1.0 / dx;

    std::vector<double> cp(grid_size - 1);
    std::vector<double> dp(grid_size);
    std::vector<double> a_imp(grid_size - 1);
    std::vector<double> b_imp(grid_size);
    std::vector<double> c_imp(grid_size - 1);
    std::vector<double> rhs(grid_size);

    // Pre-calc spatial beta
    std::vector<double> spatial_beta(grid_size);
    for (int i = 0; i < grid_size; ++i) {
        double x_curr = x_min + i * dx;
        spatial_beta[i] = (-lambda * x_curr) * drift_mult;
    }

    double diff_imp_main = 1.0 + 2.0 * C_diff;
    double diff_imp_off  = -C_diff;
    double diff_exp_main = 1.0 - 2.0 * C_diff;
    double diff_exp_off  = C_diff;
    double sum_p = 1.0 / dx;

    // --- Time Loop ---
    for (int k = 1; k < n_steps; ++k) {
        double vk = drifts[k];
        double temporal_beta = vk * drift_mult;

        double* a_ptr = a_imp.data();
        double* b_ptr = b_imp.data();
        double* c_ptr = c_imp.data();
        double* rhs_ptr = rhs.data();
        const double* p_ptr = p->data();

        std::fill(b_imp.begin(), b_imp.end(), diff_imp_main);

        for (int i = 0; i < grid_size - 1; ++i) {
            double beta_curr = temporal_beta + spatial_beta[i];
            double beta_next = temporal_beta + spatial_beta[i+1];
            c_ptr[i] = diff_imp_off + beta_next;
            a_ptr[i] = diff_imp_off - beta_curr;
        }
        b_ptr[0] = 1.0; c_ptr[0] = -1.0;
        b_ptr[grid_size - 1] = 1.0; a_ptr[grid_size - 2] = 0.0;

        rhs_ptr[0] = 0.0;
        rhs_ptr[grid_size - 1] = 0.0;
        
        for (int i = 1; i < grid_size - 1; ++i) {
            double beta_prev = temporal_beta + spatial_beta[i-1];
            double beta_next = temporal_beta + spatial_beta[i+1];
            rhs_ptr[i] = (diff_exp_off + beta_prev) * p_ptr[i-1] + 
                          diff_exp_main * p_ptr[i] + 
                         (diff_exp_off - beta_next) * p_ptr[i+1];
        }

        solve_tdma(grid_size, a_ptr, b_ptr, c_ptr, rhs, *p_new, cp, dp);

        double sum_p_new = 0.0;
        double* p_new_ptr = p_new->data();
        for(int i=0; i<grid_size; ++i) sum_p_new += p_new_ptr[i];

        double prob_absorbed = (sum_p - sum_p_new) * dx;
        if (prob_absorbed < 0) prob_absorbed = 0.0;
        likelihood[k] = prob_absorbed / dt;

        sum_p = sum_p_new;
        std::swap(p, p_new);
    }

    // --- FINAL CALCULATION: Survival Probability ---
    // Sum the probability mass remaining on the grid
    double final_mass = 0.0;
    for(double val : *p) final_mass += val;
    *survival_out = final_mass * dx;
}