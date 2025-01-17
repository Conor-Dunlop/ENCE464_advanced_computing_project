#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <chrono>
#include <iostream>

#include "poisson_mt.h"

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
    int iterations = 100;
    int n = 7;
    int threads = 4;
    float delta = 1;
    int x = -1;
    int y = -1;
    int z = -1;
    double amplitude = 1.0;
    int runs = 1;
    int opt;

    // parse the command line arguments
    while ((opt = getopt(argc, argv, "h:n:i:x:y:z:a:t:d:r:")) != -1)
    {
        switch (opt)
        {
        case 'h':
            printf("Usage: poisson [-n size] [-x source x-poisition] [-y source y-position] [-z source z-position] [-a source amplitude] [-i iterations] [-t threads] [-d] (for debug mode)\n");
            return EXIT_SUCCESS;
        case 'n':
            n = atoi(optarg);
            break;
        case 'i':
            iterations = atoi(optarg);
            break;
        case 'x':
            x = atoi(optarg);
            break;
        case 'y':
            y = atoi(optarg);
            break;
        case 'z':
            z = atoi(optarg);
            break;
        case 'a':
            amplitude = atof(optarg);
            break;
        case 't':
            threads = atoi(optarg);
            break;
        case 'd':
            debug = true;
            break;
        case 'r':
            runs = atoi(optarg);
            break;
        default:
            fprintf(stderr, "Usage: poisson [-n size] [-x source x-poisition] [-y source y-position] [-z source z-position] [-a source amplitude]  [-i iterations] [-t threads] [-d] (for debug mode)\n");
            exit(EXIT_FAILURE);
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

    // Default x,y, z
    if (x < 0 || x > n - 1)
        x = n / 2;
    if (y < 0 || y > n - 1)
        y = n / 2;
    if (z < 0 || z > n - 1)
        z = n / 2;

    source[(z * n + x) * n + y] = amplitude;

    double* result = nullptr;

    double total_time = 0;
    if (debug)
    {
        printf("Starting solver with:\n"
            "n = %i\n"
            "iterations = %i\n"
            "threads = %i\n"
            "delta = %f\n",
            n, iterations, threads, delta);
    }

    for (int i = 0; i < runs; i++)
    {
        std::chrono::time_point time_start = std::chrono::high_resolution_clock::now();
        result = poisson_mixed_multithread(n, source, iterations, threads, delta);
        std::chrono::time_point time_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> time_diff = time_end - time_start;
        total_time += time_diff.count();

        if (i < (runs-1)) {
            free(result);
        }
    }

    if (debug) {
        std::cout << "Time taken: " << total_time / (double)runs << "\n";
    }

    print_slice(n, result, 2);

    free(source);
    free(result);

    return EXIT_SUCCESS;
}

