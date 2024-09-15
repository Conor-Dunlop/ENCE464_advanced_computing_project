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


double* poisson_mixed_r2(const int n, const double* source, const int iterations, const double delta)
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

    auto accessor = [n](double* src, int i, int j, int k) {
            return src[((k+1) * (n+2) + (j+1)) * (n+2) + (i+1)];
        };

    auto cdiff_sum = [accessor](double* src, int i, int j, int k) {
        double res = accessor(src, i - 1, j, k) + accessor(src, i + 1, j, k) + accessor(src, i, j - 1, k) + accessor(src, i, j + 1, k) + accessor(src, i, j, k - 1) + accessor(src, i, j, k + 1);
        return res;
        };

    // this scares me
    for (int iter = 0; iter < iterations; iter++)
    {
        build_boundary(curr, n);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    double res = (1.0 / 6.0) * (cdiff_sum(curr, i, j, k) - (delta * delta * source[(k * n + j) * n + i]));
                    next[((k + 1) * (n + 2) + (j + 1)) * (n + 2) + (i + 1)] = res;
                }
            }
        }

        if (debug) {
            std::cout << "Iteration: " << iter << '\n';
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