/**
 * return: 0 if success, 1 otherwise
 */
int parse_args(int argc, char **argv, int *neurons_x, int *neurons_y,
               int *data_dim, double *max_epochs,
               double *initial_learning_rate,
               double *learning_rate_lambda, double *initial_radius,
               double *radius_lambda, char **input_path, char **output_path,
               int *seed);


void allocate_matrix(double **vectors, int num_vectors, int dims);


void free_matrix(double *vectors, int num_vectors, int dims);


double get_cell(double *vectors, int dims, int row, int col);


void set_cell(double *vectors, int dims, int row, int col, double val);


double *get_vec(double *vectors, int dims, int row);


/**
 * Input data must be a file whose first line is the number of vectors,
 * second line is the dimensionality,
 * all the following lines are comma-separated vector values
 */
int read_input_data(char *input_path, double **input_vectors, int *num_vectors, int *dims);


void write_output_data(char *output_path, double *output_vectors, int num_vectors, int dims);


/**
 * return: random number in range [lower, upper]
 */
long rand_in_range(long lower, long upper);


double euclidean_distance(double *v1, double *v2, int dims);
