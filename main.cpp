#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef _WIN32
    #include "getopt.h"
#else
    #include <unistd.h>
#endif

#define PRINT_THRESHOLD 15

bool debug = true;

static void build_boundary(double* buf, int32_t n) {
    const int32_t ln_s = n + 2;
    const int32_t pl_s = (n + 2) * (n + 2);
    
    for (int32_t i = 1; i <= n; i++) {
        for (int32_t j = 1; j <= n; j++) {
            buf[ln_s * (i + 1) + (j + 1) + pl_s] = -1.0; // Bottom
            buf[ln_s * (i + 1) + (j + 1) + pl_s * (n + 1)] = 1.0; // Top
        }
    }
}

double* poisson_mixed_r2(const int n, double* const source, const int iterations, const double delta) {
    if (debug) {
        printf("Starting solver with:\n"
               "n = %i\n"
               "iterations = %i\n"
               "delta = %f\n", n, iterations, delta);
    }

    size_t buf_size = (n + 2) * (n + 2) * (n + 2);
    double* curr = (double*)calloc(buf_size, sizeof(double));
    double* next = (double*)calloc(buf_size, sizeof(double));

    if (!curr || !next) {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    const double delta_sq = delta * delta;
    const double one_sixth = 1.0 / 6.0;

    for (int iter = 0; iter < iterations; iter++) {
        build_boundary(curr, n);
        for (int k = 1; k <= n; k++) {
            for (int j = 1; j <= n; j++) {
                double* row_b = curr + ((k * (n + 2) + j) * (n + 2));
                double* row_f = curr + ((k * (n + 2) + j) * (n + 2) + 2);
                double* src_row = source + ((k - 1) * n + (j - 1) * n);

                for (int i = 0; i < n; i++) {
                    double res = row_b[i] + row_f[i] + 
                                 curr[((k + 1) * (n + 2) + j) * (n + 2) + i + 1] +
                                 curr[((k - 1) * (n + 2) + j) * (n + 2) + i + 1] +
                                 curr[(k * (n + 2) + j + 1) * (n + 2) + i + 1] +
                                 curr[(k * (n + 2) + j - 1) * (n + 2) + i + 1];

                    next[(k * (n + 2) + j) * (n + 2) + i + 1] = one_sixth * (res - (delta_sq * src_row[i]));
                }
            }
        }

        double* temp = curr;
        curr = next;
        next = temp;
    }

    free(next);
    return curr;
}

void print_slice(int n, double* data) {
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            printf("%0.5f ", data[((n / 2) * (n + 2) + y + 1) * (n + 2) + (x + 1)]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    int iterations = 10, n = 5;
    double delta = 1;
    double amplitude = 1.0;
    int x = -1, y = -1, z = -1;

    int opt;
    while ((opt = getopt(argc, argv, "h:n:i:x:y:z:a:d:")) != -1) {
        switch (opt) {
            case 'h':
                printf("Usage: poisson [-n size] [-x source x-poisition] [-y source y-position] [-z source z-position] [-a source amplitude] [-i iterations] [-d] (for debug mode)\n");
                return EXIT_SUCCESS;
            case 'n': n = atoi(optarg); break;
            case 'i': iterations = atoi(optarg); break;
            case 'x': x = atoi(optarg); break;
            case 'y': y = atoi(optarg); break;
            case 'z': z = atoi(optarg); break;
            case 'a': amplitude = atof(optarg); break;
            case 'd': debug = true; break;
            default: fprintf(stderr, "Usage: poisson [-n size] [-x source x-poisition] [-y source y-position] [-z source z-position] [-a source amplitude] [-i iterations] [-d] (for debug mode)\n"); exit(EXIT_FAILURE);
        }
    }

    if (n % 2 == 0) {
        fprintf(stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    double* source = (double*)calloc(n * n * n, sizeof(double));
    if (!source) {
        fprintf(stderr, "Error: failed to allocate source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    if (x < 0) x = n / 2;
    if (y < 0) y = n / 2;
    if (z < 0) z = n / 2;

    source[(z * n + y) * n + x] = amplitude;

    double* result_r2 = poisson_mixed_r2(n, source, iterations, delta);

    if (debug) {
        printf("Finished solving.\n");
    }

    if (n <= PRINT_THRESHOLD) {
        print_slice(n, result_r2);
    }

    free(source);
    free(result_r2);
    return EXIT_SUCCESS;
}
