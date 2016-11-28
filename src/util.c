#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"


unsigned long flop_counter = 0;


enum FlagType {
    NONE,
    NEURONS_X,
    NEURONS_Y,
    DATA_DIM,
    INPUT_PATH,
    OUTPUT_PATH,
    MAX_EPOCHS,
    INITIAL_A,
    INITIAL_R,
    LAMBDA_A,
    LAMBDA_R,
    PRINT_FLOPS
};


int parse_int(char *str, int *target_int) {
    char *str_end;

    long dummy = strtol(str, &str_end, 10);

    if (str != str_end) {
        *target_int = (int) dummy;
        return 0;
    }

    printf("Could not convert %s to int\n", str);

    return 1;
}


int parse_float(char *str, double *target_float) {
    char *str_end;

    double dummy = strtod(str, &str_end);

    if (str != str_end) {
        *target_float = (double) dummy;
        return 0;
    }

    printf("Could not convert %s to float\n", str);

    return 1;
}

int get_flag(char *option) {
    if (strcmp(option, "-x") == 0) {
        return NEURONS_X;
    } else if (strcmp(option, "-y") == 0) {
        return NEURONS_Y;
    } else if (strcmp(option, "-m") == 0) {
        return DATA_DIM;
    } else if (strcmp(option, "-i") == 0) {
        return INPUT_PATH;
    } else if (strcmp(option, "-o") == 0) {
        return OUTPUT_PATH;
    } else if (strcmp(option, "-e") == 0) {
        return MAX_EPOCHS;
    } else if (strcmp(option, "-a0") == 0) {
        return INITIAL_A;
    } else if (strcmp(option, "-al") == 0) {
        return LAMBDA_A;
    } else if (strcmp(option, "-r0") == 0) {
        return INITIAL_R;
    } else if (strcmp(option, "-rl") == 0) {
        return LAMBDA_R;
    } else if (strcmp(option, "-p") == 0) {
        return PRINT_FLOPS;
    } else {
        return NONE;
    }
}


int parse_args(int argc,
                char **argv,
                int *neurons_x,
                int *neurons_y,
                int *data_dim,
                double *max_epochs,
                double *initial_learning_rate,
                double *learning_rate_lambda,
                double *initial_radius,
                double *radius_lambda,
                char **input_path,
                char **output_path,
                int *print_flops) {

    char *cmd_part;

    int last_flag = NONE, success = 0;

    for (int i = 0; i < argc; ++i) {
        cmd_part = argv[i];

        if (cmd_part[0] == '-') {
            last_flag = get_flag(cmd_part);

            if (last_flag == PRINT_FLOPS) {
                *print_flops = 1;

                last_flag = NONE;
            }

            continue;
        }

        // This command part is an option value!

        switch (last_flag) {
            case NONE: // No flag was set at all
                {
                    break;
                }
            case NEURONS_X:
                {
                    success = parse_int(cmd_part, neurons_x);

                    break;
                }
            case NEURONS_Y:
                {
                    success = parse_int(cmd_part, neurons_y);

                    break;
                }
            case DATA_DIM:
                {
                    success = parse_int(cmd_part, data_dim);

                    break;
                }
            case INPUT_PATH:
                {
                    *input_path = cmd_part;

                    break;
                }
            case OUTPUT_PATH:
                {
                    *output_path = cmd_part;

                    break;
                }
            case MAX_EPOCHS:
                {
                    success = parse_float(cmd_part, max_epochs);

                    break;
                }
            case INITIAL_A:
                {
                    success = parse_float(cmd_part, initial_learning_rate);

                    break;
                }
            case INITIAL_R:
                {
                    success = parse_float(cmd_part, initial_radius);

                    break;
                }
            case LAMBDA_R:
                {
                    success = parse_float(cmd_part, radius_lambda);

                    break;
                }
            case LAMBDA_A:
                {
                    success = parse_float(cmd_part, learning_rate_lambda);

                    break;
                }
            default:
                {
                    success = 1;    // This should never run.
                }
        }

        if (success != 0) {
            return success;
        }

        last_flag = NONE;

    }

    return success;

}


void allocate_matrix(double **vectors, int num_vectors, int dims) {
    *vectors = (double *) calloc(num_vectors * dims, sizeof(double));
}


void free_matrix(double *vectors, int num_vectors, int dims) {
    free(vectors);
}


double get_cell(double *vectors, int dims, int row, int col) {
    return vectors[row * dims + col];
}


void set_cell(double *vectors, int dims, int row, int col, double val) {
    vectors[row * dims + col] = val;
}


double *get_vec(double *vectors, int dims, int row) {
    return vectors + row * dims;
}


/**
 * Return success indicator
 */
int parse_vector(char *line, double *vector, int dims) {
    char *remaining = line;
    double dummy;

    // While we have things to parse
    for (int i = 0; i < dims; ++i) {
        dummy = strtod(line, &remaining);

        // No more things to parse but there should be
        if (remaining == line) {
            return 1;
        }

        vector[i] = dummy;

        // There must be either a comma or a newline
        if (remaining == '\0') {
            return 1;
        }

        line = remaining + 1;
    }

    return 0;
}


int read_input_data(char *input_path, double **input_vectors, int *num_vectors, int *dims) {
    int success = 0;

    FILE *input_file = fopen(input_path, "r");

    int MAX_LINE_SIZE = 1024, vector_index = 0;

    char line[MAX_LINE_SIZE];

    for (int line_number = 0; fgets(line, MAX_LINE_SIZE, input_file) != NULL; ++line_number) {
        if (line_number == 0) {
            // Read number of vectors
            success = parse_int(line, num_vectors);

            if (success != 0) {
                fclose(input_file);
                return success;
            }
        } else if (line_number == 1) {
            // Read vector dimensions
            success = parse_int(line, dims);

            if (success != 0) {
                fclose(input_file);
                return success;
            }

            // Allocate matrix
            allocate_matrix(input_vectors, *num_vectors, *dims);
        } else {
            vector_index = line_number - 2;

            // Read vector
            success = parse_vector(line, get_vec(*input_vectors, *dims, vector_index), *dims);

            if (success != 0) {
                free_matrix(*input_vectors, *num_vectors, *dims);

                fclose(input_file);
                return success;
            }
        }
    }

    fclose(input_file);
    return success;
}


void write_output_data(char *output_path, double *output_vectors, int num_vectors, int dims) {
    FILE *output_file = fopen(output_path, "w");

    if (dims > 0) {
        fprintf(output_file, "%d\n", num_vectors);
        fprintf(output_file, "%d\n", dims);
        for (int vector_num = 0; vector_num < num_vectors; ++vector_num) {
            fprintf(output_file, "%f", get_cell(output_vectors, dims, vector_num, 0));

            for (int dim = 1; dim < dims; ++dim) {
                fprintf(output_file, ",%f", get_cell(output_vectors, dims, vector_num, dim));
            }

            fprintf(output_file, "\n");
        }
    }

    fclose(output_file);
}


/**
 * Taken from an answer by Ryan Reich (http://stackoverflow.com/a/6852396)
 */
long random_at_most(long max) {
  unsigned long
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    num_bins = (unsigned long) max + 1,
    num_rand = (unsigned long) RAND_MAX + 1,
    bin_size = num_rand / num_bins,
    defect   = num_rand % num_bins;

  long x;
  do {
   x = random();
  }
  // This is carefully written not to overflow
  while (num_rand - defect <= (unsigned long)x);

  // Truncated division is intentional
  return x/bin_size;
}


long rand_in_range(long lower, long upper) {
    return lower + random_at_most(upper - lower);
}


double euclidean_distance(double *v1, double *v2, int dims) {
    double sum = 0.0;

    for (int i = 0; i < dims; ++i) {
        flop_counter += 4;
        sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }

    ++flop_counter;
    return sqrt(sum);
}
