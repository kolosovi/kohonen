#include <stdio.h>

#include "util.h"


int ROWS = 1000, COLS = 2000;


float_type ** test_allocate_matrix() {
    float_type **vectors;

    allocate_matrix(&vectors, ROWS, COLS);

    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            vectors[i][j] = 123456789.0;
        }
    }

    return vectors;
}


void test_free_matrix(float_type **vectors, int num_vectors, int dims) {
    free_matrix(vectors, num_vectors, dims);
}


void test_read_input_data(char *filename) {
    int num_vectors = 0, dims = 0;
    float_type **vectors;
    float_type *vector;

    read_input_data(filename, &vectors, &num_vectors, &dims);

    for (int i = 0; i < num_vectors; ++i) {
        vector = vectors[i];

        for (int j = 0; j < dims; ++j) {
            printf("%f ", vector[j]);
        }

        printf("\n");
    }
}

void test_rand_in_range() {
    for (int i = 0; i < 10; ++i) {
        printf("%f ", (float_type) rand_in_range(0, 255));
    }
}

int main(int argc, char **argv) {
    float_type **matrix;

    matrix = test_allocate_matrix();

    test_free_matrix(matrix, ROWS, COLS);

    test_read_input_data("test_data/input_matrix.txt");

    test_rand_in_range();
}
