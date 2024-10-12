#include "calc.h"
#include <immintrin.h>

void compute_func(const int k, double* curr, double* next, double* const source, const int n, const double delta) {

    const int sz = n + 2;

    const double delta_sq_six = delta * delta * (1.0 / 6.0);
    __m512d one_sixth = _mm512_set1_pd(1.0 / 6.0);
    __m512d dsq_sixth = _mm512_set1_pd(delta_sq_six);

    const int iter_step = 8;
    int final_row_sz = n%iter_step;
    int num_iter = n/iter_step;

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