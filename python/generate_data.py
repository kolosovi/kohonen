import argparse
import sys

import numpy as np


CLUSTER_CENTERS = [
    np.array([255.0, 0.0, 0.0]),    # Red
    np.array([0.0, 255.0, 0.0]),    # Green
    np.array([0.0, 0.0, 255.0]),    # Blue
]


def generate_data(cluster_size, standard_deviation, output):
    data = []

    for center in CLUSTER_CENTERS:
        for i in xrange(cluster_size):
            noise = np.random.normal(0, standard_deviation, 3)
            data.append(np.clip(
                center + noise,
                0.0,
                255.0
            ))

    output.write("{data_size}\n{data_dim}\n".format(
        data_size=(cluster_size * 3),
        data_dim=3
    ))

    for item in data:
        output.write("{0},{1},{2}\n".format(
            item[0], item[1], item[2]
        ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate input data"
    )

    parser.add_argument("-c", "--cluster-size",
                        metavar="NUMBER",
                        type=int,
                        help="Cluster size",
                        required=True
                       )
    parser.add_argument("-o", "--output",
                        metavar="PATH",
                        type=str,
                        help="Output path"
                       )
    parser.add_argument("-s", "--standard-deviation",
                        metavar="FLOAT",
                        type=float,
                        help="Standard deviation of points in a cluster",
                        default=20.0
                       )
    args = parser.parse_args()

    if args.output is not None:
        output = open(args.output, "w")
    else:
        output = sys.stdout

    generate_data(args.cluster_size, args.standard_deviation, output)
