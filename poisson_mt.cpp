#include "poisson_mt.h"

static void build_boundary(double* buf, int32_t n) {
    const int32_t ln_s = (n + 2);
    const int32_t pl_s = (n + 2) * (n + 2);
    const int32_t top_plane_idx = pl_s * (n + 1);
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++)
        {
            buf[ln_s * (i + 1) + (j + 1) + pl_s] = -1.0;
            buf[ln_s * (i + 1) + (j + 1) + pl_s * (n)] = 1.0;
            buf[pl_s * (i + 1) + ln_s * (j + 1)] = buf[pl_s * (i + 1) + ln_s * (j + 1) + 2];
            buf[pl_s * (i + 1) + ln_s * (j + 1) + (n + 1)] = buf[pl_s * (i + 1) + ln_s * (j + 1) + (n - 1)];
            buf[pl_s * (i + 1) + (j + 1)] = buf[pl_s * (i + 1) + (j + 1) + ln_s * 2];
            buf[pl_s * (i + 1) + (j + 1) + ln_s * (n + 1)] = buf[pl_s * (i + 1) + (j + 1) + ln_s * (n - 1)];
        }
    }
}

static void compute_func(const int k, double* curr, double* next, double* const source, const int n, const double delta) {

    const int sz = n + 2;

    const double delta_sq_six = delta * delta * (1.0 / 6.0);
    __m512d one_sixth = _mm512_set1_pd(1.0 / 6.0);
    __m512d dsq_sixth = _mm512_set1_pd(delta_sq_six);

    const int iter_step = 8;
    int final_row_sz = n%iter_step;
    int num_iter = n/iter_step;

    __mmask8 defmask = (unsigned char)0xff;
    __mmask8 finalmask = (unsigned char)(0xff >> (n - final_row_sz));
    
    __mmask8 mask = _load_mask8(&defmask);

    for (int j = 1; j < n + 1; j++)
    {
        double* row_b = curr + ((k * sz + j) * sz - 0);
        double* row_f = curr + ((k * sz + j) * sz + 2);
        double* row_up = curr + (((k + 1) * sz + j) * sz + 1);
        double* row_low = curr + (((k - 1) * sz + j) * sz + 1);
        double* row_l = curr + ((k * sz + (j + 1)) * sz + 1);
        double* row_r = curr + ((k * sz + (j - 1)) * sz + 1);
        double* src_row = source + (((k - 1) * n + (j - 1)) * n - 0);
        double* dest_row = next + ((k * sz + j) * sz + 1);

        

        for (int i = 0; i < n; i += iter_step)
        {
            __m512d a = _mm512_load_pd(&row_b[i]);
            __m512d b = _mm512_load_pd(&row_f[i]);
            __m512d res_a = _mm512_add_pd(a, b);

            __m512d c = _mm512_load_pd(&row_up[i]);
            __m512d d = _mm512_load_pd(&row_low[i]);
            __m512d res_b = _mm512_add_pd(c, d);

            __m512d e = _mm512_load_pd(&row_l[i]);
            __m512d f = _mm512_load_pd(&row_r[i]);
            __m512d res_c = _mm512_add_pd(e, f);
            
            __m512d res_mid = _mm512_add_pd(res_a, res_b);
            __m512d res = _mm512_add_pd(res_mid, res_c);

            __m512d src_val_int = _mm512_load_pd(&src_row[i]);
            __m512d src_val = _mm512_mul_pd(src_val_int, dsq_sixth);
            __m512d final = _mm512_fmsub_pd(res, one_sixth, src_val);
            _mm512_store_pd(&dest_row[i], final);
        }
    }
}

double* poisson_mixed_multithread(const int n, double* const source, const int iterations, const int threads, const double delta)
{
    if (debug)
    {
        printf("Starting solver with:\n"
            "n = %i\n"
            "iterations = %i\n"
            "threads = %i\n"
            "delta = %f\n",
            n, iterations, threads, delta);
    }

    // Allocate some buffers to calculate the solution in
    size_t buf_size = (n + 2) * (n + 2) * (n + 2);

    double* curr = (double*)calloc(buf_size, sizeof(double));
    double* next = (double*)calloc(buf_size, sizeof(double));

    // Ensure we haven't run out of memory
    if (curr == nullptr || next == nullptr)
    {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    dp::thread_pool pool(threads);
    
    for (int iter = 0; iter < iterations; iter++)
    {
        build_boundary(curr, n);

        for (int k = 1; k < n + 1; k++)
        {
            pool.enqueue_detach(compute_func, k, curr, next, source, n, delta);
        }
        
        pool.wait_for_tasks();

        // swap current and next
        double* curr_temp = curr;
        curr = next;
        next = curr_temp;
    }


    // Free one of the buffers and return the correct answer in the other.
    // The caller is now responsible for free'ing the returned pointer.
    free(next);

    if (debug)
    {
        printf("Finished solving.\n");
    }

    return curr;
}