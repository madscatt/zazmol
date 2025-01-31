#include <stdio.h>

void matrix_multiply_c(int r1, int c1, int c2, float *first, float *second, float *result) {
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            result[i * c2 + j] = 0.0;
            for (int k = 0; k < c1; ++k) {
                result[i * c2 + j] += first[i * c1 + k] * second[k * c2 + j];
            }
        }
    }
}
