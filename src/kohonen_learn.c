#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <mpi.h>

// Custom headers
#include "util.h"


static int num_tasks, rank;
static int neuron_num_start = 0;       // Included in the range
static int neuron_num_end = 0;         // Not included in the range
static int process_num_neurons = 0;    // # of neurons given to this process

static double start_time = 0.0, end_time = 0.0;

enum ProcessType {
    MAIN_PROCESS
};


/**
 * Initialize neuron weight vectors to random values
 * This should be done in a centralized fashion, and the neuron weights should
 * be broadcasted to all processes. Maybe the same should be done on every
 * iteration of the learning process. Probably not.
 */
void kohonen_init_weights(double *weight_vectors,
                     int data_dim,
                     int num_neurons) {
    for (int i = 0; i < num_neurons; ++i) {
        for (int j = 0; j < data_dim; ++j) {
            set_cell(weight_vectors, data_dim, i, j, (double) rand_in_range(0, 255));
        }
    }
}


/**
 * What is the output of learning? There's no output probably. We just have
 * a double** matrix that represents neuron weights. This matrix gets
 * changed on every iteration.
 */
int kohonen_learn_iter(double *input_vectors,
                  int input_size,
                  int data_dim,
                  double *weight_vectors,
                  int num_neurons,
                  int neurons_x,
                  int neurons_y,
                  double epoch,
                  double initial_learning_rate,
                  double learning_rate_lambda,
                  double initial_radius,
                  double radius_lambda) {

    double *input_vector;

    // Min dist between input vector and neurons for each process
    double bmu_dists[num_tasks];
    double min_distance = -1.0;
    double radius = 0.0, learning_rate = 0.0;
    double lattice_dists[process_num_neurons];

    bmu_dists[rank] = -1.0;

    int bmu_numbers[num_tasks];

    int bmu_number = 0, bmu_x = 0, bmu_y = 0, vector_number = 0;
    int global_neuron_number = 0;
    int neuron_x = 0, neuron_y = 0;
    int mpi_success = 0;
    int min_distance_rank = 0;

    // 1. Let main process pick a random input vector
    if (rank == MAIN_PROCESS) {
        vector_number = (int) rand_in_range(0, input_size - 1);
    }

    mpi_success = MPI_Bcast((void *)(&vector_number), 1, MPI_INT, MAIN_PROCESS, MPI_COMM_WORLD);

    if (mpi_success != 0) {
        return mpi_success;
    }

    // Each process gets the same input vector
    input_vector = get_vec(input_vectors, data_dim, vector_number);

    // 2. Calculate distances from neuron weight vectors to the input vector
    // and find the neuron with minimum distance (the best-matching unit, BMU)
    for (int neuron_number = 0; neuron_number < process_num_neurons; ++neuron_number) {
        double *weight_vector = get_vec(weight_vectors, data_dim, neuron_number);

        double distance = euclidean_distance(
            input_vector,
            weight_vector,
            data_dim);

        if (bmu_dists[rank] < 0.0 || distance < bmu_dists[rank]) {
            bmu_dists[rank] = distance;
            bmu_numbers[rank] = neuron_num_start + neuron_number;
        }
    }

    // Mpi gather distances and local BMU numbers and THEN find the minimum.
    mpi_success = MPI_Gather(
        (void *) (bmu_dists + rank),
        1,
        MPI_DOUBLE,
        (void *) bmu_dists,
        1,
        MPI_DOUBLE,
        MAIN_PROCESS,
        MPI_COMM_WORLD);

    if (mpi_success != 0) {
        return mpi_success;
    }

    mpi_success = MPI_Gather(
        (void *) (bmu_numbers + rank),
        1,
        MPI_INT,
        (void *) bmu_numbers,
        1,
        MPI_INT,
        MAIN_PROCESS,
        MPI_COMM_WORLD);

    if (mpi_success != 0) {
        return mpi_success;
    }

    // Main process finds BMU number
    if (rank == MAIN_PROCESS) {
        for (int i = 0; i < num_tasks; ++i) {
            if (min_distance < 0.0 || bmu_dists[i] < min_distance) {
                min_distance = bmu_dists[i];
                min_distance_rank = i;
            }
        }

        bmu_number = bmu_numbers[min_distance_rank];
    }

    // Main process broadcasts BMU number to all processes
    mpi_success = MPI_Bcast(
        (void *)(&bmu_number),
        1,
        MPI_INT,
        MAIN_PROCESS,
        MPI_COMM_WORLD);

    bmu_x = bmu_number / neurons_y;
    bmu_y = bmu_number % neurons_y;

    // 3. Calculate neighborhood radius and learning rate
    radius = initial_radius * exp(-(epoch / radius_lambda));
    learning_rate = initial_learning_rate *
                    exp(-(epoch / learning_rate_lambda));

    // 4. Calculate distances from BMU to other neurons on the lattice
    // neuron_x = neuron_num_start / neurons_y;
    // neuron_y = neuron_num_start % neurons_y;
    for (int neuron_number = 0; neuron_number < process_num_neurons; ++neuron_number) {
        global_neuron_number = neuron_num_start + neuron_number;
        neuron_x = global_neuron_number / neurons_y;
        neuron_y = global_neuron_number % neurons_y;

        if (global_neuron_number == bmu_number) {
            lattice_dists[neuron_number] = 0.0;
        } else {
            lattice_dists[neuron_number] = (double) sqrt(
                    (neuron_x - bmu_x) * (neuron_x - bmu_x) +
                    (neuron_y - bmu_y) * (neuron_y - bmu_y));
        }

        // neuron_y = (neuron_y + 1) % neurons_y;
        // if (neuron_y == 0) {
        //     ++neuron_x;
        // }
    }

    // 5. Re-calculate weight vectors
    for (int neuron_number = 0; neuron_number < process_num_neurons; ++neuron_number) {
        double d = lattice_dists[neuron_number];

        if (d > radius) {
            continue;
        }

        for (int dim = 0; dim < data_dim; ++dim) {
            double w = get_cell(weight_vectors, data_dim, neuron_number, dim);
            double x = input_vector[dim];

            w = w + exp(-d / (2.0 * pow(radius, 2))) * learning_rate * (x - w);

            set_cell(weight_vectors, data_dim, neuron_number, dim, w);
        }
    }

    return 0;
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
 *   ./kohonen_learn -i test_input_data.txt -x 40 -y 40 -m 3 -e 1000 -a0 0.1 -al 100.0 -r0 40.0 -rl 100.0 -o neuron_weights.txt works nice
 */
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int success = 0;

    success = MPI_Barrier(MPI_COMM_WORLD);

    if (success != 0) {
        MPI_Finalize();
        return 1;
    }

    start_time = MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc <= 1) {
        if (rank == MAIN_PROCESS) {
            print_help(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int num_vectors = 0;
    double *input_vectors;
    double *local_neuron_weights;
    double *neuron_weights = NULL;  // Used only in main process

    int neurons_x = -1, neurons_y = -1, data_dim = -1, seed = -1, num_neurons = 0;

    double max_epochs = -1.0, initial_learning_rate = -1.0, learning_rate_lambda = -1.0, initial_radius = -1.0, radius_lambda = -1.0;

    char *input_path = "FOOBAR", *output_path = "FOOBAR";

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

    srandom(rank * time(NULL));

    num_neurons = neurons_x * neurons_y;

    success = read_input_data(input_path, &input_vectors, &num_vectors, &data_dim);

    if (success != 0) {
        MPI_Finalize();
        return 1;
    }

    // Number of neurons each task handles
    int neuron_range_lengths[num_tasks];

    // Neuron range start * data dimensions
    int neuron_displacements[num_tasks];

    int neuron_range_length = num_neurons / num_tasks;
    int neuron_remainder = num_neurons % num_tasks;

    neuron_num_start = rank * neuron_range_length;

    if (rank < neuron_remainder) {
        neuron_num_start += rank;
        neuron_num_end = neuron_num_start + neuron_range_length + 1;
    } else {
        neuron_num_start += neuron_remainder;
        neuron_num_end = neuron_num_start + neuron_range_length;
    }

    if (rank == (num_tasks - 1)) {
        neuron_num_end = num_neurons;
    }

    neuron_displacements[rank] = neuron_num_start * data_dim;

    process_num_neurons = (neuron_num_end - neuron_num_start);
    neuron_range_lengths[rank] = process_num_neurons * data_dim;

    // Gather neuron range lengths so we can gather neuron weights later
    success =  MPI_Gather(
                    (void *) (neuron_range_lengths + rank),
                    1, MPI_INT,
                    neuron_range_lengths,
                    1, MPI_INT,
                    MAIN_PROCESS,
                    MPI_COMM_WORLD);

    if (success != 0) {
        MPI_Finalize();
        return 1;
    }

    // Gather neuron displacements so we can gather neuron weights later
    success =  MPI_Gather(
                    (void *) (neuron_displacements + rank),
                    1, MPI_INT,
                    neuron_displacements,
                    1, MPI_INT,
                    MAIN_PROCESS,
                    MPI_COMM_WORLD);

    if (success != 0) {
        MPI_Finalize();
        return 1;
    }

    // Initialize neuron weights. We broadcast them to all processes, but
    // this should not be necessary because each process operates on its
    // own range of neurons anyway
    allocate_matrix(&local_neuron_weights, process_num_neurons, data_dim);

    // Let each process initialize the weights of its neurons
    kohonen_init_weights(local_neuron_weights, data_dim, process_num_neurons);

    for (double epoch = 0; epoch < max_epochs; epoch += 1.0) {
        success = kohonen_learn_iter(input_vectors,
                                     num_vectors,
                                     data_dim,
                                     local_neuron_weights,
                                     num_neurons,
                                     neurons_x,
                                     neurons_y,
                                     epoch,
                                     initial_learning_rate,
                                     learning_rate_lambda,
                                     initial_radius,
                                     radius_lambda);

        if (success != 0) {
            MPI_Finalize();
            return success;
        }
    }

    // Gather neuron weights from all processes to the main process. This is pretty tough......
    // This requires us to first gather data range lengths, then use gatherv to gather neuron weights
    // Then the main process writes neuron weights to output path
    if (rank == MAIN_PROCESS) {
        allocate_matrix(&neuron_weights, num_neurons, data_dim);
    }
    success = MPI_Gatherv(
        (void *)(local_neuron_weights),
        neuron_range_lengths[rank],
        MPI_DOUBLE,
        (void *)neuron_weights,
        neuron_range_lengths,
        neuron_displacements,
        MPI_DOUBLE,
        MAIN_PROCESS,
        MPI_COMM_WORLD);

    if (success != 0) {
        MPI_Finalize();
        return 1;
    }

    if (rank == MAIN_PROCESS) {
        write_output_data(output_path, neuron_weights, num_neurons, data_dim);
    }

    success = MPI_Barrier(MPI_COMM_WORLD);

    if (success != 0) {
        MPI_Finalize();
        return 1;
    }

    end_time = MPI_Wtime();

    if (rank == MAIN_PROCESS) {
        printf(" %f\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
