from casacore.tables import table
import numpy as np
from argparse import ArgumentParser


def euclidean_distance(pos1, pos2):
    """Euclidean distance"""
    return np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2)


def parse_args():
    """
    Parse input arguments
    """

    parser = ArgumentParser(description='Longest baseline length')
    parser.add_argument('ms', help='MS')
    parser.add_argument('--min_dist', default=1500, type=float, help='Minimal distance in km')
    return parser.parse_args()


def main():
    """
    Main script
    """

    args = parse_args()

    # Open the ANTENNA table of your Measurement Set
    ms_antenna = table(f'{args.ms}/ANTENNA')

    # Get antenna positions and names
    antenna_positions = ms_antenna.getcol('POSITION')
    antenna_names = ms_antenna.getcol('NAME')

    # Initialize variables to track the longest baseline
    longest_baseline = 0
    longest_baseline_antennas = None

    # Calculate all possible baselines and track the longest
    print(f"Baselines beyond {args.min_dist} km:")
    n_antennas = len(antenna_positions)
    for i in range(n_antennas):
        for j in range(i + 1, n_antennas):
            distance = euclidean_distance(antenna_positions[i], antenna_positions[j])
            if distance > longest_baseline:
                longest_baseline = distance
                longest_baseline_antennas = (antenna_names[i], antenna_names[j])
            if round(distance/1000, 2)>args.min_dist:
                print(f"{round(distance/1000, 2)} km between {antenna_names[i]} and {antenna_names[j]}")

    # Output the longest baseline and the antennas forming it
    print(f"\nLongest Baseline: {round(longest_baseline/1000, 2)} km")
    print(f"Formed by antennas: {longest_baseline_antennas[0]} and {longest_baseline_antennas[1]}")


if __name__ == '__main__':
    main()
