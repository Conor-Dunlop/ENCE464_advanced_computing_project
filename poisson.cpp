#include "poisson.h"

double* poisson_mixed(const int n, const double* source, const int iterations, const double delta)
{
    if (debug)
    {
        printf("Starting solver with:\n"
            "n = %i\n"
            "iterations = %i\n"
            "delta = %f\n",
            n, iterations, delta);
    }

    // Allocate some buffers to calculate the solution in
    size_t buffer_size = n * n * n * sizeof(double);

    double* curr = (double*)calloc(n * n * n, sizeof(double));
    double* next = (double*)calloc(n * n * n, sizeof(double));

    // Ensure we haven't run out of memory
    if (curr == nullptr || next == nullptr)
    {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    auto accessor = [n](double* src, int i, int j, int k) {
        // enforce boundary conditions and extrapolate points for central difference
        auto sgn = [](int32_t t) { return (t < 0) ? -1 : (t > 0); };

        bool i_valid = i >= 0 && i < n;
        bool j_valid = j >= 0 && j < n;
        bool k_valid = k >= 1 && k < n - 1;

        // branching go brrrr
        if (i_valid && j_valid && k_valid) {
            return src[(k * n + j) * n + i];
        }

        else if (!i_valid) {
            return src[(k * n + j) * n + (i - 2*sgn(i))]; // neumann boundary
        }

        else if (!j_valid) {
            return src[(k * n + (j - 2 * sgn(j))) * n + i];
        }

        else if (!k_valid) {
            return 1.0 * (double)sgn(k - 1); // dirichlet boundary
        }

        else {
            return 0.0; // should never happen
        }

        };

    auto cdiff_sum = [accessor, curr](int i, int j, int k) {
        double res = accessor(curr, i - 1, j, k) + accessor(curr, i + 1, j, k) + accessor(curr, i, j - 1, k) + accessor(curr, i, j + 1, k) + accessor(curr, i, j, k - 1) + accessor(curr, i, j, k + 1);
        return res;
        };

    // this scares me
    for (int iter = 0; iter < iterations; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    double res = (1.0 / 6.0) * (cdiff_sum(i, j, k) - (delta * delta * source[(k * n + j) * n + i]));
                    next[(k * n + j) * n + i] = res;
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