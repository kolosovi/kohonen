#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Custom headers
#include "util.h"


/**
 * Initialize neuron weight vectors to random values
 */
void kohonen_init_weights(float_type **weight_vectors,
                     int data_dim,
                     int num_neurons) {
    for (int i = 0; i < num_neurons; ++i) {
        for (int j = 0; j < data_dim; ++j) {
            weight_vectors[i][j] = (float_type) rand_in_range(0, 255);
        }
    }
}


/**
 * What is the output of learning? There's no output probably. We just have
 * a double** matrix that represents neuron weights. This matrix gets
 * changed on every iteration.
 */
void kohonen_learn_iter(float_type **input_vectors,
                   int input_size,
                   int data_dim,
                   float_type **weight_vectors,
                   int num_neurons,
                   int neurons_x,
                   int neurons_y,
                   float_type epoch,
                   float_type initial_learning_rate,
                   float_type learning_rate_lambda,
                   float_type initial_radius,
                   float_type radius_lambda) {
    // printf("Epoch %f\n", epoch); fflush(stdout);
    float_type *input_vector;
    float_type min_distance = -1.0, distance = 0.0;
    float_type radius = 0.0, learning_rate = 0.0;
    float_type lattice_dists[num_neurons];

    int bmu_number = 0, bmu_x = 0, bmu_y = 0;
    int vector_number = (int) rand_in_range(0, input_size - 1);

    // printf("Vector index %d\n", vector_number); fflush(stdout);

    // 1. Pick random input vector
    input_vector = input_vectors[vector_number];

    // 2. Calculate distances from neuron weight vectors to the input vector
    // and find the neuron with minimum distance (the best-matching unit, BMU)
    for (int x = 0; x < neurons_x; ++x) {
        for (int y = 0; y < neurons_y; ++y) {
            int neuron_number = (x * neurons_y) + y;

            distance = euclidean_distance(
                input_vector,
                weight_vectors[neuron_number],
                data_dim);

            if (min_distance < 0.0 || distance < min_distance) {
                min_distance = distance;
                bmu_number = neuron_number;
                bmu_x = x;
                bmu_y = y;
            }
        }
    }

    // 3. Calculate neighborhood radius and learning rate
    radius = initial_radius * exp(-(epoch / radius_lambda));
    learning_rate = initial_learning_rate *
                    exp(-(epoch / learning_rate_lambda));

    printf("Learning rate %f\n", learning_rate); fflush(stdout);
    printf("Radius %f\n", radius); fflush(stdout);

    // 4. Calculate distances from BMU to other neurons on the lattice
    for (int x = 0; x < neurons_x; ++x) {
        for (int y = 0; y < neurons_y; ++y) {
            int neuron_number = (x * neurons_y) + y;

            if (neuron_number == bmu_number) {
                lattice_dists[neuron_number] = 0.0;
            } else {
                lattice_dists[neuron_number] = (float_type) sqrt(
                        (x - bmu_x) * (x - bmu_x) +
                        (y - bmu_y) * (y - bmu_y));
            }
        }
    }

    // 5. Re-calculate weight vectors
    for (int i = 0; i < num_neurons; ++i) {
        float_type d = lattice_dists[i];

        if (d > radius) {
            continue;
        }

        for (int dim = 0; dim < data_dim; ++dim) {
            float_type w = weight_vectors[i][dim];
            float_type x = input_vector[dim];

            // printf("%d distance: %f\n", i, d); fflush(stdout);

            // float_type modifier = exp(-d / (2.0 * pow(radius, 2))) * learning_rate * (x - w);

            // printf("%d term: %f\n", i, modifier); fflush(stdout);

            w = w + exp(-d / (2.0 * pow(radius, 2))) * learning_rate * (x - w);

            weight_vectors[i][dim] = w;
        }
    }
}


void print_help(char *command_name) {
    char *help = "Usage: %s options, where options are:\n"
    "    -i - Path to input data file (required)\n"
    "    -x - Number of neurons by x\n"
    "    -y - Number of neurons by y\n"
    "    -m - Data dimensionality\n"
    "    -e - Max number of epochs\n"
    "    -a0 - Initial learning rate\n"
    "    -al - Learning rate lambda parameter\n"
    "    -r0 - Initial radius\n"
    "    -rl - Radius lambda parameter\n"
    "    -s - PRNG seed\n";

    printf(help, command_name);
}


/**
 *   -x -- Number of neurons by x
 *   -y -- Number of neurons by y
 *   -m -- Data dimensionality
 *   -i -- Path to input data file
 *   -e -- Max number of epochs
 *   -a0 -- Initial learning rate
 *   -al -- Learning rate lambda parameter
 *   -r0 -- Initial radius
 *   -rl -- Radius lambda parameter
 *   -s -- PRNG seed (optional)
 */
int main(int argc, char **argv) {
    if (argc <= 1) {
        print_help(argv[0]);
        return 1;
    }

    int success = 0;

    int num_vectors = 0;
    float_type **input_vectors;
    float_type **neuron_weights;

    int neurons_x = -1, neurons_y = -1, data_dim = -1, seed = -1, num_neurons = 0;

    float_type max_epochs = -1.0, initial_learning_rate = -1.0, learning_rate_lambda = -1.0, initial_radius = -1.0, radius_lambda = -1.0;

    char *input_path = "FOOBAR", *output_path = "FOOBAR";

    srandom(time(NULL));

    parse_args(argc, argv,
                &neurons_x,
                &neurons_y,
                &data_dim,
                &max_epochs,
                &initial_learning_rate,
                &learning_rate_lambda,
                &initial_radius,
                &radius_lambda,
                &input_path,
                &output_path,
                &seed);

    num_neurons = neurons_x * neurons_y;

    printf(
        "X: %d, Y: %d, dim: %d, seed: %d\n",
        neurons_x, neurons_y, data_dim, seed
    );

    printf(
        "Epochs: %f, initial a: %f, lambda a: %f, initial radius: %f, lambda radius: %f\n",
        max_epochs, initial_learning_rate, learning_rate_lambda, initial_radius, radius_lambda
    );

    printf(
        "Input path: %s\n",
        input_path
    );

    success = read_input_data(input_path, &input_vectors, &num_vectors, &data_dim);

    if (success != 0) {
        exit(1);
    }

    // Initialize neuron weights
    allocate_matrix(&neuron_weights, num_neurons, data_dim);
    kohonen_init_weights(neuron_weights, data_dim, num_neurons);

    for (float_type epoch = 0; epoch < max_epochs; epoch += 1.0) {
        kohonen_learn_iter(input_vectors,
                           num_vectors,
                           data_dim,
                           neuron_weights,
                           num_neurons,
                           neurons_x,
                           neurons_y,
                           epoch,
                           initial_learning_rate,
                           learning_rate_lambda,
                           initial_radius,
                           radius_lambda);
    }

    write_output_data(output_path, neuron_weights, num_neurons, data_dim);
}
