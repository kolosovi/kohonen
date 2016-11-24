import argparse

from matplotlib import pyplot as plt
import numpy as np

def clip(v):
    return min(max(0, int(v)), 255)


def get_weights(input_path, neurons_x, neurons_y):
    weight_list = []

    n_weights = np.zeros((neurons_x, neurons_y, 3), dtype=np.uint8)

    with open(input_path, "r") as nfile:
        for i, line in enumerate(nfile):
            if i > 1:
                weights = map(float, line.strip().split(","))
                weight_list.append(weights)

    for i in xrange(neurons_x):
        for j in xrange(neurons_y):
            num = i * neurons_x + j
            w_arr = weight_list[num]
            n_weights[i][j][0] = clip(w_arr[0])
            n_weights[i][j][1] = clip(w_arr[1])
            n_weights[i][j][2] = clip(w_arr[2])

    return n_weights


def plot_weights(input_path, neurons_x, neurons_y):
    n_weights = get_weights(input_path, neurons_x, neurons_y)
    
    plt.imsave("image.png", n_weights, vmin=0.0, vmax=255.0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot neuron weights"
    )

    parser.add_argument("-i", "--input",
                        metavar="PATH",
                        type=str,
                        help="Input path"
                       )
    parser.add_argument("-x", "--neurons_x",
                        metavar="NUMBER",
                        type=int,
                        help="Number of neurons by x"
                       )
    parser.add_argument("-y", "--neurons_y",
                        metavar="NUMBER",
                        type=int,
                        help="Number of neurons by y"
                       )
    args = parser.parse_args()

    plot_weights(args.input, args.neurons_y, args.neurons_x)
