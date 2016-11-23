#include <stdio.h>
#include <stdlib.h>

#include "util.h"


int ROWS = 1000, COLS = 2000;


double * test_allocate_matrix() {
    double *vectors;

    allocate_matrix(&vectors, ROWS, COLS);

    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            set_cell(vectors, COLS, i, j, 123456789.0);
        }
    }

    return vectors;
}


void test_read_input_data(char *filename) {
    int num_vectors = 0, dims = 0;
    double *vectors;
    double *vector;

    read_input_data(filename, &vectors, &num_vectors, &dims);

    for (int i = 0; i < num_vectors; ++i) {
        vector = get_vec(vectors, dims, i);

        for (int j = 0; j < dims; ++j) {
            printf("%f ", vector[j]);
        }

        printf("\n");
    }
}

void test_rand_in_range() {
    for (int i = 0; i < 10; ++i) {
        printf("%f ", (double) rand_in_range(0, 255));
    }
}

int main(int argc, char **argv) {
    double *matrix;

    matrix = test_allocate_matrix();

    test_read_input_data("test_data/input_matrix.txt");

    test_rand_in_range();

    free(matrix);
}
