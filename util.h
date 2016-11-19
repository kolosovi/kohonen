#include "base.h"


/**
 * return: 0 if success, 1 otherwise
 */
int parse_args(int argc, char **argv, int *neurons_x, int *neurons_y,
               int *data_dim, float_type *max_epochs,
               float_type *initial_learning_rate,
               float_type *learning_rate_lambda, float_type *initial_radius,
               float_type *radius_lambda, char **input_path, char **output_path,
               int *seed);


void allocate_matrix(float_type ***vectors, int num_vectors, int dims);


void free_matrix(float_type **vectors, int num_vectors, int dims);


/**
 * Input data must be a file whose first line is the number of vectors,
 * second line is the dimensionality,
 * all the following lines are comma-separated vector values
 */
int read_input_data(char *input_path, float_type ***input_vectors, int *num_vectors, int *dims);


void write_output_data(char *output_path, float_type **output_vectors, int num_vectors, int dims);


/**
 * return: random number in range [lower, upper]
 */
long rand_in_range(long lower, long upper);


float_type euclidean_distance(float_type *v1, float_type *v2, int dims);
