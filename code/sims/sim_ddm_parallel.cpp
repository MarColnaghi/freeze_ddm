#include "mex.h"
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include <omp.h> // Header for OpenMP

// --- 1. XorShift128+ (Same as before) ---
struct FastRNG {
    uint64_t s[2];
    
    // Constructor now vital for seeding per thread
    FastRNG(uint64_t seed) {
        // SplitMix64 to robustly initialize state from a 64-bit seed
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

// --- 2. Marsaglia Polar Method (Same as before) ---
struct FastGaussian {
    FastRNG rng;
    double cache;
    bool has_cache;
    const double SCALE; 

    FastGaussian(uint64_t seed) : rng(seed), cache(0.0), has_cache(false), SCALE(2.0 / 9007199254740992.0) {}

    inline double next() {
        if (has_cache) {
            has_cache = false;
            return cache;
        }
        double u, v, s;
        do {
            u = (rng.next() >> 11) * SCALE - 1.0;
            v = (rng.next() >> 11) * SCALE - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);

        double mul = std::sqrt(-2.0 * std::log(s) / s);
        cache = v * mul;
        has_cache = true;
        return u * mul;
    }
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 3) mexErrMsgIdAndTxt("sim_ddm_par:invalidNumInputs", "Inputs: drifts, params, n_trials, [seed]");

    double* drifts = mxGetPr(prhs[0]);
    size_t n_time_steps = mxGetNumberOfElements(prhs[0]);

    double* params = mxGetPr(prhs[1]);
    double dt = params[0];
    double sigma_sq = params[1];
    double x0 = params[2];
    double theta = params[3];
    double lambda = params[4];
    int n_trials = (int)mxGetScalar(prhs[2]);

    // Optional seed input
    uint64_t base_seed = 12345;
    if (nrhs > 3) {
        base_seed = (uint64_t)mxGetScalar(prhs[3]);
    }

    plhs[0] = mxCreateDoubleMatrix(n_trials, 1, mxREAL);
    double* fpt_results = mxGetPr(plhs[0]);

    double sigma = std::sqrt(sigma_sq);
    double noise_scale = sigma * std::sqrt(dt);
    double decay = 1.0 - lambda * dt;
    
    // Pre-calc drifts (Shared read-only memory is fine for threads)
    std::vector<double> drifts_dt(n_time_steps);
    for(size_t i=0; i<n_time_steps; ++i) {
        drifts_dt[i] = drifts[i] * dt;
    }
    const double* d_dt_ptr = drifts_dt.data();

    // --- OPENMP PARALLEL SECTION ---
    // This directive tells the compiler to spawn threads and split the loop
    #pragma omp parallel
    {
        // 1. Setup Thread-Local RNG
        // Each thread gets a unique seed based on its ID
        int thread_id = omp_get_thread_num();
        FastGaussian local_gen(base_seed + thread_id * 9999); 

        // 2. The Parallel Loop
        // "schedule(static)" divides the work into equal chunks (fastest for equal workloads)
        #pragma omp for schedule(static)
        for (int i = 0; i < n_trials; ++i) {
            double x = x0;
            double finish_time = std::numeric_limits<double>::quiet_NaN();

            // Inner loop (Time steps) - Runs locally in each thread
            for (int k = 0; k < n_time_steps; ++k) {
                x = x * decay + d_dt_ptr[k] + noise_scale * local_gen.next();

                if (x >= theta) {
                    finish_time = (k + 1) * dt; 
                    break; 
                }
            }
            fpt_results[i] = finish_time;
        }
    }
}