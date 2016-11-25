"""
1st dimension is number of epochs, 2d dimension is number of processes, 3d
dimension is execution time.
"""
import argparse
import os
import re

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


NEURON_FLOP = {
    4000: 2018078725,
    2000: 1005706195,
    1000: 503746486,
    8000: 4017177166,
    6000: 3026981284,
    9000: 4533140599,
    10000: 5035779622,
    5000: 2520140860,
    7000: 3528285415,
    3000: 1514202850
}


EPOCH_FLOP = {
    10000: 656136343,
    20000: 1315127714,
    40000: 2635371310,
    50000: 3299596985,
    30000: 1977046869
}

# Parameter types
NEURON = "neuron"
EPOCH = "epoch"

PARAMETER_LABELS = {
    NEURON: "# of neurons",
    EPOCH: "# of iterations"
}

# Plot types
EXECUTION_TIME = "time"
PARALLELIZATION_EFFICIENCY = "parallelization efficiency"
PERFORMANCE = "performance"
EFFICIENCY = "efficiency"

VALUE_LABELS = {
    EXECUTION_TIME: "Execution time",
    PARALLELIZATION_EFFICIENCY: "Parallelization efficiency",
    PERFORMANCE: "Performance",
    EFFICIENCY: "Efficiency"
}

# Other constants
LOMONOSOV_PEAK_PERFORMANCE = 1.7 * (1e15)


name_re = re.compile(r"stdout_(?P<num_epochs>\d+)_(epochs|neurons)_(?P<num_procs>\d+)_processes.txt")


def file_list(path):
    return [os.path.join(path, name) for name in os.listdir(path) if "stdout" in name]


def show_plot(path, parameter, plot_type):
    filenames = file_list(path)
    X, Y, Z = get_x_y_z(filenames, parameter, plot_type)

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(PARAMETER_LABELS[parameter])
    ax.set_ylabel("# of processes")
    ax.set_xlabel(VALUE_LABELS[plot_type])
    ax.plot_surface(X, Y, Z, shade=False, rstride=1, cstride=1, cmap=cm.coolwarm)

    plt.show()

def get_x_y_z(input_paths, parameter, plot_type):
    xy_to_z = {}

    for input_path in input_paths:
        match = re.search(name_re, input_path)

        if not match:
            continue

        num_epochs = int(match.groupdict()["num_epochs"])
        num_procs = int(match.groupdict()["num_procs"])

        with open(input_path, "r") as process_stdout:
            execution_time = float(process_stdout.read().strip())

        # This is another row to the matrix
        xy_to_z[(num_epochs, num_procs)] = execution_time

    epochs = np.asarray(sorted(set(k[0] for k in xy_to_z.iterkeys())))
    procs = np.asarray(sorted(set(k[1] for k in xy_to_z.iterkeys())))

    X, Y = np.meshgrid(epochs, procs)
    Z = np.zeros(X.shape)

    proc_dims = X.shape[0]
    epoch_dims = X.shape[1]

    for epoch_i in xrange(epoch_dims):
        for proc_i in xrange(proc_dims):
            num_epochs = X[0][epoch_i]
            num_procs = Y[proc_i][0]

            Z[proc_i, epoch_i] = xy_to_z[(num_epochs, num_procs)]

    # Z holds execution time now.

    if plot_type == PARALLELIZATION_EFFICIENCY:
        Z = (Z[0] / Z) / Y
    elif plot_type == PERFORMANCE:
        # Divide FLOP by time
        if parameter == NEURON:
            flop = NEURON_FLOP
        elif parameter == EPOCH:
            flop = EPOCH_FLOP

        for epoch_i in xrange(epoch_dims):
            epoch_flops = flop[X[0][epoch_i]]

            for proc_i in xrange(proc_dims):
                print "Divide value for {0} {1}, {2} procs by {3}".format(
                    X[0][epoch_i],
                    parameter,
                    Y[proc_i][0],
                    epoch_flops
                )
                Z[proc_i, epoch_i] /= epoch_flops
    elif plot_type == EFFICIENCY:
        # Efficiency is performance divided by peak performance
        pass

    return X, Y, Z


def dump_plot(input_paths, output_path):
    X, Y, Z = get_x_y_z(input_paths)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_wireframe(X, Y, Z)

    plt.savefig(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Construct execution time matrix"
    )

    parser.add_argument("-i", "--input",
                        metavar="PATH",
                        type=str,
                        nargs="+",
                        help="Input paths"
                       )
    parser.add_argument("-o", "--output",
                        metavar="PATH",
                        type=str,
                        help="Output path"
                       )
    args = parser.parse_args()

    dump_plot(args.input, args.output)

