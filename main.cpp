#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <chrono>
#include <iostream>

#include "poisson.h"
#include "poisson_r2.h"
#include "poisson_multithread.h"

#define PRINT_THRESHOLD 15

/**
 * poisson.c
 * Implementation of a Poisson solver with Dirichlet boundary conditions.
 *
 * This template handles the basic program launch, argument parsing, and memory
 * allocation required to implement the solver *at its most basic level*. You
 * will likely need to allocate more memory, add threading support, account for
 * cache locality, etc...
 *
 * BUILDING:
 * gcc -o poisson poisson.c -lpthread
 *
 * [note: linking pthread isn't strictly needed until you add your
 *        multithreading code]
 *
 * TODO:
 * 1 - Read through this example, understand what it does and what it gives you
 *     to work with.
 * 2 - Implement the basic algorithm and get a correct output.
 * 3 - Add a timer to track how long your execution takes.
 * 4 - Profile your solution and identify weaknesses.
 * 5 - Improve it!
 * 6 - Remember that this is now *your* code and *you* should modify it however
 *     needed to solve the assignment.
 *
 * See the lab notes for a guide on profiling and an introduction to
 * multithreading (see also threads.c which is reference by the lab notes).
 */


 // Global flag
 // Set to true when operating in debug mode to enable verbose logging
bool debug = false;

void print_slice(int n, double* data, int ex_size = 0) {

    uint32_t val_add = ex_size / 2;

    for (int x = 0; x < n; ++x)
    {
        for (int y = 0; y < n; ++y)
        {
            // printf("%0.5f ", result[((n / 2) * n + y) * n + x]);
            printf("%0.5f ", data[(((n + ex_size) / 2) * (n + ex_size) + (y + val_add)) * (n + ex_size) + (x+val_add)]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv)
{
    // Default settings for solver
    int iterations = 300;
    int n = 101;
    int threads = 1;
    double delta = 1.0;

#ifdef _DEBUG || DEBUG
    debug = true;
#endif // DEBUG

    // parse the command line arguments
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            printf("Usage: poisson [-n size] [-i iterations] [-t threads] [--debug]\n");
            return EXIT_SUCCESS;
        }

        if (strcmp(argv[i], "-n") == 0)
        {
            if (i == argc - 1)
            {
                fprintf(stderr, "Error: expected size after -n!\n");
                return EXIT_FAILURE;
            }

            n = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "-i") == 0)
        {
            if (i == argc - 1)
            {
                fprintf(stderr, "Error: expected iterations after -i!\n");
                return EXIT_FAILURE;
            }

            iterations = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "-t") == 0)
        {
            if (i == argc - 1)
            {
                fprintf(stderr, "Error: expected threads after -t!\n");
                return EXIT_FAILURE;
            }

            threads = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "--debug") == 0)
        {
            debug = true;
        }
    }

    // Ensure we have an odd sized cube
    if (n % 2 == 0)
    {
        fprintf(stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    // Create a source term with a single point in the centre
    double* source = (double*)calloc(n * n * n, sizeof(double));
    if (source == NULL)
    {
        fprintf(stderr, "Error: failed to allocated source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    source[(n * n * n) / 2] = 1;

    ///////////////////////////////////////////////////////////////////////////////

    std::chrono::time_point time_start_r1 = std::chrono::high_resolution_clock::now();

    // Calculate the resulting field with Dirichlet conditions
    double* result_r1 = poisson_mixed(n, source, iterations, delta);

    std::chrono::time_point time_end_r1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_diff_r1 = time_end_r1 - time_start_r1;
    std::cout << "Time taken for R1: " << time_diff_r1 << '\n';
    if (n <= PRINT_THRESHOLD) {
        std::cout << "Result:" << '\n';
        print_slice(n, result_r1);
    }
    std::cout << '\n';

    ///////////////////////////////////////////////////////////////////////////////

    std::chrono::time_point time_start_r2 = std::chrono::high_resolution_clock::now();

    // Calculate the resulting field with Dirichlet conditions
    double* result_r2 = poisson_mixed_r2(n, source, iterations, delta);

    std::chrono::time_point time_end_r2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_diff_r2 = time_end_r2 - time_start_r2;
    std::cout << "Time taken for R2: " << time_diff_r2 << "\n";
    if (n <= PRINT_THRESHOLD) {
        std::cout << "Result:" << '\n';
        print_slice(n, result_r2, 2);
    }
    std::cout << '\n';


    free(source);
    free(result_r1);
    free(result_r2);

    return EXIT_SUCCESS;
}

