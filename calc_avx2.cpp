#include "calc.h"
#include <immintrin.h>

#include "calc.h"
#include <immintrin.h>

inline __m256i get_mask(const uint8_t mask){
    // https://stackoverflow.com/questions/21622212/how-to-perform-the-inverse-of-mm256-movemask-epi8-vpmovmskb
    __m256i vmask(_mm256_set1_epi64x(mask));
    const __m256i admask(_mm256_setr_epi64x(0x1, 0x2, 0x4, 0x8));
    vmask = _mm256_and_si256(vmask, admask);
    const __m256i bit_mask(_mm256_setr_epi64x(0xfffffffffffffffe, 0xfffffffffffffffd, 0xfffffffffffffffb, 0xfffffffffffffff7));
    vmask = _mm256_or_si256(vmask, bit_mask);
    return _mm256_cmpeq_epi64(vmask, _mm256_set1_epi64x(-1));
}

void compute_func(const int k, double* curr, double* next, double* const source, const int n, const double delta) {

    const int sz = n + 2;

    const double delta_sq_six = delta * delta * (1.0 / 6.0);
    __m256d one_sixth = _mm256_set1_pd(1.0 / 6.0);
    __m256d dsq_sixth = _mm256_set1_pd(delta_sq_six);

    const int iter_step = 4; // Eight doubles in one AVX512 register

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
        
        size_t i = 0;
        for (; i+iter_step <= n; i+=iter_step)
        {
            __m256d a = _mm256_load_pd(&row_b[i]);
            __m256d b = _mm256_load_pd(&row_f[i]);
            __m256d res_a = _mm256_add_pd(a, b);

            __m256d c = _mm256_load_pd(&row_up[i]);
            __m256d d = _mm256_load_pd(&row_low[i]);
            __m256d res_b = _mm256_add_pd(c, d);

            __m256d e = _mm256_load_pd(&row_l[i]);
            __m256d f = _mm256_load_pd(&row_r[i]);
            __m256d res_c = _mm256_add_pd(e, f);
            
            __m256d res_mid = _mm256_add_pd(res_a, res_b);
            __m256d res = _mm256_add_pd(res_mid, res_c);

            __m256d src_val_int = _mm256_load_pd(&src_row[i]);
            __m256d src_val = _mm256_mul_pd(src_val_int, dsq_sixth);
            __m256d final = _mm256_fmsub_pd(res, one_sixth, src_val);
            _mm256_store_pd(&dest_row[i], final);
        }

        if (i < n) // The remaining number is less than 8
        {
            __m256i final_mask = get_mask((uint8_t)((1<<(n-i))-1));

            __m256d a = _mm256_maskload_pd(&row_b[i], final_mask);
            __m256d b = _mm256_maskload_pd(&row_f[i], final_mask);
            __m256d res_a = _mm256_add_pd(a, b);

            __m256d c = _mm256_maskload_pd(&row_up[i], final_mask);
            __m256d d = _mm256_maskload_pd(&row_low[i], final_mask);
            __m256d res_b = _mm256_add_pd(c, d);

            __m256d e = _mm256_maskload_pd(&row_l[i], final_mask);
            __m256d f = _mm256_maskload_pd(&row_r[i], final_mask);
            __m256d res_c = _mm256_add_pd(e, f);
            
            __m256d res_mid = _mm256_add_pd(res_a, res_b);
            __m256d res = _mm256_add_pd(res_mid, res_c);

            __m256d src_val_int = _mm256_maskload_pd(&src_row[i], final_mask);
            __m256d src_val = _mm256_mul_pd(src_val_int, dsq_sixth);
            __m256d final = _mm256_fmsub_pd(res, one_sixth, src_val);
            _mm256_maskstore_pd(&dest_row[i], final_mask, final);
        }
    }
}