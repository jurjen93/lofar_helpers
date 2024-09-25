import pandas as pd
import numpy as np
from math import radians, cos, sin, sqrt, atan2
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.cm as cm
import matplotlib.pyplot as plt


def angular_distance(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """
    Function to calculate angular distance between two points with RA/DEC

    :param ra1_deg: RA first source
    :param dec1_deg: DEC first source
    :param ra2_deg: RA second source
    :param dec2_deg: DEC second source

    :return: distance in degrees

    """
    # Convert degrees to radians
    ra1 = radians(ra1_deg)
    dec1 = radians(dec1_deg)
    ra2 = radians(ra2_deg)
    dec2 = radians(dec2_deg)

    # Haversine-like formula for angular distance in degrees
    delta_ra = ra2 - ra1
    delta_dec = dec2 - dec1

    a = sin(delta_dec / 2) ** 2 + cos(dec1) * cos(dec2) * sin(delta_ra / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    # Convert back to degrees
    distance_deg = np.degrees(c)
    return distance_deg


def remove_close_pairs(csv_file, dist_offset=0.1):
    """
    Remove close pairs in CSV and sort by scalarphasediff score

    Args:
        csv_file: CSV file with sources RA/DEC and spd_score
        dist_offset: min offset in degrees

    Returns:
        cleaned df
    """

    # Load the CSV file
    df = pd.read_csv(csv_file).sort_values('spd_score')

    # Assuming the CSV file has columns 'RA' and 'DEC' in degrees
    ra = df['RA'].values
    dec = df['DEC'].values

    # List to store pairs and their distances
    pairs = []

    # Calculate pairwise distances
    for i in range(len(ra)):
        for j in range(i + 1, len(ra)):
            dist = angular_distance(ra[i], dec[i], ra[j], dec[j])
            if dist < dist_offset:
                pairs.append((i, j, dist))

    # Convert to a DataFrame to display
    duplicate_df = pd.DataFrame(pairs, columns=['Index 1', 'Index 2', 'Distance (degrees)'])

    for idx in duplicate_df['Index 2'].values[::-1]:
        df.drop(idx, inplace=True)

    return df


def plot_voronoi(csv_path, ra_column='RA', dec_column='DEC', score_column='spd_score'):
    """
    Plots a Voronoi diagram where the points (RA/DEC) are colored based on the values in the score_column.

    Parameters:
        csv_path (str): Path to the CSV file.
        ra_column (str): Column name for RA values.
        dec_column (str): Column name for DEC values.
        score_column (str): Column name for the score values used for coloring the points.
    """
    # Load the CSV data
    data = pd.read_csv(csv_path).sort_values(score_column)

    # Extract RA, DEC, and spd_scores columns
    ra = data[ra_column].values
    dec = data[dec_column].values
    spd_scores = data[score_column].values

    # Combine RA and DEC into coordinate pairs
    points = np.column_stack((ra, dec))

    # Compute Voronoi tessellation
    vor = Voronoi(points)

    # Create a plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the Voronoi diagram without filling the regions
    voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors='black', line_width=2, line_alpha=0.6)

    # Normalize spd_scores for coloring
    norm = plt.Normalize(vmin=0, vmax=2)
    cmap = cm.viridis

    # Scatter plot the points, colored by spd_scores
    sc = ax.scatter(ra[1:], dec[1:], c=spd_scores[1:], cmap=cmap, norm=norm, edgecolor='black', s=200, zorder=2)

    # Highlight the first point with a red star
    ax.scatter(ra[0], dec[0], color='red', marker='*', s=600, zorder=3, edgecolor='black')


    ax.tick_params(labelsize=16)

    # Add a colorbar
    # Add a colorbar and set font sizes
    cbar = plt.colorbar(sc)
    cbar.set_label('$\hat{\sigma}_{c}$ (rad)', fontsize=16)  # Set the font size for the colorbar label
    cbar.ax.tick_params(labelsize=16)  # Set the font size for colorbar ticks
    plt.gca().invert_xaxis()

    # Set labels and title
    ax.set_xlabel('RA (degrees)', fontsize=16)
    ax.set_ylabel('DEC (degrees)', fontsize=16)

    plt.tight_layout()
    plt.savefig('voronoi_calibrators.png', dpi=150)