#include "calc.h"

void compute_func(const int k, double* curr, double* next, double* const source, const int n, const double delta) {

    const int sz = n + 2;
    const double one_sixth = (1.0 / 6.0);
    const double delta_sq = delta * delta;

    
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

        for (int i = 0; i < n; i++)
        {
            double res = (
                row_b[i] +
                row_f[i] +
                row_up[i] +
                row_low[i] +
                row_l[i] +
                row_r[i]);

            double src_val = one_sixth * (res - (delta_sq * src_row[i]));
            dest_row[i] = src_val;
        }
    }
}
