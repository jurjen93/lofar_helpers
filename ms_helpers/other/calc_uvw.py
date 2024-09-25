from casacore.tables import table
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, GCRS
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

# Ensure the latest Earth Orientation Parameters are used
iers.conf.auto_download = True

def recalculate_uvw(ms_path):
    """
    Recalculate UVW coordinates for a measurement set.

    Parameters:
    - ms_path (str): Path to the measurement set.

    Returns:
    - None: The function updates the UVW coordinates in the MS in-place.
    """

    # Open the measurement set in writable mode
    with table(ms_path, readonly=False, ack=False) as ms:
        # Read necessary columns from the MAIN table
        time_col = ms.getcol('TIME')          # Observation times in seconds since MJD epoch
        antenna1_col = ms.getcol('ANTENNA1')  # Indices of first antennas
        antenna2_col = ms.getcol('ANTENNA2')  # Indices of second antennas

        # Get antenna positions from the ANTENNA table
        with table(ms_path + '::ANTENNA', ack=False) as ant_table:
            ant_positions = ant_table.getcol('POSITION')  # Shape: (n_antennas, 3)
            n_antennas = ant_table.nrows()

        # Get source direction from the FIELD table
        with table(ms_path + '::FIELD', ack=False) as field_table:
            phase_dir = field_table.getcol('PHASE_DIR')  # Shape: (nfields, 1, 2)
            ra = phase_dir[0, 0, 0]  # Right Ascension in radians
            dec = phase_dir[0, 0, 1]  # Declination in radians

        # Convert times to astropy Time objects
        times = Time(time_col / 86400.0, format='mjd', scale='utc')

        # Convert antenna positions to EarthLocation objects
        ant_locations = EarthLocation(x=ant_positions[:, 0] * u.m,
                                      y=ant_positions[:, 1] * u.m,
                                      z=ant_positions[:, 2] * u.m)

        # Initialize arrays to hold GCRS positions
        ant_gcrs = []

        # Transform antenna positions from ITRF to GCRS at each time
        # Since the antennas are fixed on Earth, their positions in GCRS change over time due to Earth's rotation.

        print("Transforming antenna positions to GCRS frame...")
        for ant_loc in ant_locations:
            itrs = ant_loc.get_itrs(obstime=times[0])
            gcrs = itrs.transform_to(GCRS(obstime=times[0]))
            ant_gcrs.append(gcrs.cartesian)

        ant_gcrs = np.array([[gcrs.x.value, gcrs.y.value, gcrs.z.value] for gcrs in ant_gcrs])  # Shape: (n_antennas, 3)

        # Compute baselines in GCRS frame
        bl_gcrs = ant_gcrs[antenna2_col] - ant_gcrs[antenna1_col]  # Shape: (nrows, 3)

        # Define the source direction in GCRS frame
        source_icrs = SkyCoord(ra=ra * u.rad, dec=dec * u.rad, frame='icrs')
        source_gcrs = source_icrs.transform_to(GCRS(obstime=times[0]))  # Use the first time for transformation

        # Get the unit vector pointing to the source in GCRS frame
        src_unit_vector = np.array([source_gcrs.cartesian.x.value,
                                    source_gcrs.cartesian.y.value,
                                    source_gcrs.cartesian.z.value])  # Shape: (3,)

        # Create orthonormal basis (u_vec, v_vec, w_vec) where w_vec is towards the source
        w_vec = src_unit_vector / np.linalg.norm(src_unit_vector)  # Normalize w_vec

        # Choose an arbitrary vector not parallel to w_vec to construct u_vec
        arbitrary_vector = np.array([0, 0, 1])
        if np.allclose(w_vec, arbitrary_vector):
            arbitrary_vector = np.array([0, 1, 0])

        # Compute u_vec = cross(arbitrary_vector, w_vec)
        u_vec = np.cross(arbitrary_vector, w_vec)
        u_vec /= np.linalg.norm(u_vec)  # Normalize u_vec

        # Compute v_vec = cross(w_vec, u_vec)
        v_vec = np.cross(w_vec, u_vec)

        # Now project the baselines onto the (u_vec, v_vec, w_vec) basis
        uvw = np.zeros_like(bl_gcrs)
        uvw[:, 0] = np.dot(bl_gcrs, u_vec)
        uvw[:, 1] = np.dot(bl_gcrs, v_vec)
        uvw[:, 2] = np.dot(bl_gcrs, w_vec)

        # Update the UVW column in the measurement set
        ms.putcol('UVW', uvw)

        print("UVW recalculation completed.")
