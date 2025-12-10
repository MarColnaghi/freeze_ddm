#include "mex.h"
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

// --- 1. XorShift128+ (Same Fast RNG) ---
struct FastRNG {
    uint64_t s[2];
    
    FastRNG(uint64_t seed) {
        uint64_t z = (seed + 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        s[0] = z ^ (z >> 31);
        z = (s[0] + 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        s[1] = z ^ (z >> 31);
    }

    inline uint64_t next() {
        uint64_t s1 = s[0];
        const uint64_t s0 = s[1];
        s[0] = s0;
        s1 ^= s1 << 23; 
        s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); 
        return s[1] + s0;
    }
};

// --- 2. Marsaglia Polar Method (No Trigonometry) ---
struct FastGaussian {
    FastRNG rng;
    double cache;
    bool has_cache;
    // Constant for converting uint64 to double [-1, 1)
    const double SCALE; 

    FastGaussian(uint64_t seed) : rng(seed), cache(0.0), has_cache(false), SCALE(2.0 / 9007199254740992.0) {}

    inline double next() {
        if (has_cache) {
            has_cache = false;
            return cache;
        }

        double u, v, s;
        // Rejection loop: typically runs 1.27 times on average (very fast)
        do {
            // Generate u, v in [-1, 1) directly from bits
            // We avoid the double conversion overhead twice by doing it here
            u = (rng.next() >> 11) * SCALE - 1.0;
            v = (rng.next() >> 11) * SCALE - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);

        // One sqrt and one log, NO sin/cos
        double mul = std::sqrt(-2.0 * std::log(s) / s);
        
        cache = v * mul;
        has_cache = true;
        return u * mul;
    }
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 3) mexErrMsgIdAndTxt("sim_ddm_ultra:invalidNumInputs", "Inputs: drifts, params, n_trials");

    double* drifts = mxGetPr(prhs[0]);
    size_t n_time_steps = mxGetNumberOfElements(prhs[0]);

    double* params = mxGetPr(prhs[1]);
    double dt = params[0];
    double sigma_sq = params[1];
    double x0 = params[2];
    double theta = params[3];
    double lambda = params[4];
    int n_trials = (int)mxGetScalar(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(n_trials, 1, mxREAL);
    double* fpt_results = mxGetPr(plhs[0]);

    double sigma = std::sqrt(sigma_sq);
    double noise_scale = sigma * std::sqrt(dt);
    
    // Pre-calculate decay factor to replace multiplication with subtraction
    // x_new = x_old + (v - lambda*x)*dt 
    // x_new = x_old * (1 - lambda*dt) + v*dt
    double decay = 1.0 - lambda * dt;
    
    // Pre-multiply drifts by dt to remove multiplication inside loop
    std::vector<double> drifts_dt(n_time_steps);
    for(size_t i=0; i<n_time_steps; ++i) {
        drifts_dt[i] = drifts[i] * dt;
    }
    const double* d_dt_ptr = drifts_dt.data();

    FastGaussian gen(12345); 

    for (int i = 0; i < n_trials; ++i) {
        double x = x0;
        double finish_time = std::numeric_limits<double>::quiet_NaN();

        for (int k = 0; k < n_time_steps; ++k) {
            // Optimization: Reduced math operations per step
            // Old: x += (drifts[k] - lambda * x) * dt + scale * rand;
            // New: x = x * decay + drift_dt[k] + scale * rand;
            
            x = x * decay + d_dt_ptr[k] + noise_scale * gen.next();

            if (x >= theta) {
                finish_time = (k + 1) * dt; 
                break; 
            }
        }
        fpt_results[i] = finish_time;
    }
}