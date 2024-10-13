#include "poisson_mt.h"

static void build_boundary(double* buf, int32_t n) {
    const int32_t ln_s = (n + 2);
    const int32_t pl_s = (n + 2) * (n + 2);
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

double* poisson_mixed_multithread(const int n, double* const source, const int iterations, const int threads, const double delta)
{
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
        printf("Finished solving\n");
    }

    return curr;
}