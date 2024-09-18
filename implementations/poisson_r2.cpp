#include "poisson_r2.h"

void build_boundary(double* buf, int32_t n) {
    const int32_t ln_s = (n + 2);
    const int32_t pl_s = (n + 2) * (n + 2);
    const int32_t top_plane_idx = pl_s * (n+1);
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

void test_build_boundary(void) {
    uint32_t n = 5;

    size_t buf_size = (n + 2) * (n + 2) * (n + 2);
    double* curr = (double*)calloc(buf_size, sizeof(double));

    int32_t val = 0.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                curr[((k + 1) * (n + 2) + (j + 1)) * (n + 2) + (i + 1)] = (double)val;
                ++val;
            }
        }
    }

    build_boundary(curr, n);
}


double* poisson_mixed_r2(const int n, double* const source, const int iterations, const double delta)
{
    if (debug)
    {
        printf("Starting solver with:\n"
            "n = %i\n"
            "iterations = %i\n"
            "delta = %f\n",
            n, iterations, delta);

        // test_build_boundary();
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

    // this scares me
    const int sz = n + 2;
    const double one_sixth = (1.0 / 6.0);
    const double delta_sq = delta * delta;
    for (int iter = 0; iter < iterations; iter++)
    {
        build_boundary(curr, n);

        for (int k = 1; k < n+1; k++)
        {
            for (int j = 1; j < n+1; j++)
            {
                double* row_b = curr + ((k * sz + j) * sz - 1);
                double* row_f = curr + ((k * sz + j) * sz + 1);
                double* row_up = curr + (((k+1) * sz + j) * sz);
                double* row_low = curr + (((k-1) * sz + j) * sz);
                double* row_l = curr + ((k * sz + (j+1)) * sz);
                double* row_r = curr + ((k * sz + (j-1)) * sz);
                double* src_row = source + (((k - 1) * n + (j - 1)) * n - 1);

                for (int i = 1; i < n+1; i++)
                {

                    double res = (
                        row_b[i] +
                        row_f[i] +
                        row_up[i] +
                        row_low[i] +
                        row_l[i] +
                        row_r[i]);

                    double src_val = (delta_sq * src_row[i]);
                    next[(k * sz + j) * sz + i] = one_sixth * (res - src_val);
                }
            }
        }

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