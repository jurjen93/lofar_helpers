import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5


def compute_julian_century(epoch):
    """
    Compute the Julian century for a given epoch relative to J2000.0.
    """
    JD = epoch.jd
    JD2000 = 2451545.0  # Julian date for J2000.0
    T = (JD - JD2000) / 36525.0
    return T


def compute_precession_angles(epoch1, epoch2):
    """
    Compute the precession angles to transform coordinates from epoch1 to epoch2.
    Based on the IAU 2000 precession model.
    """
    T1 = compute_julian_century(epoch1)
    T2 = compute_julian_century(epoch2)
    dT = T2 - T1

    zeta_A = (2306.2181 + 1.39656 * T1 - 0.000139 * T1 ** 2) * dT \
             + (0.30188 - 0.000344 * T1) * dT ** 2 + 0.017998 * dT ** 3
    theta_A = (2004.3109 - 0.85330 * T1 - 0.000217 * T1 ** 2) * dT \
              - (0.42665 + 0.000217 * T1) * dT ** 2 - 0.041833 * dT ** 3
    z_A = (2306.2181 + 1.39656 * T1 - 0.000139 * T1 ** 2) * dT \
          + (1.09468 + 0.000066 * T1) * dT ** 2 + 0.018203 * dT ** 3

    # Convert to radians
    zeta_A = np.deg2rad(zeta_A / 3600)
    theta_A = np.deg2rad(theta_A / 3600)
    z_A = np.deg2rad(z_A / 3600)

    return zeta_A, theta_A, z_A


def rotation_matrix(angle, axis):
    """
    Create a rotation matrix for a given angle and axis.
    """
    if axis == 'x':
        return np.array([[1, 0, 0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    elif axis == 'y':
        return np.array([[np.cos(angle), 0, np.sin(angle)],
                         [0, 1, 0],
                         [-np.sin(angle), 0, np.cos(angle)]])
    elif axis == 'z':
        return np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0],
                         [0, 0, 1]])


def compute_precession_matrix(epoch1, epoch2):
    """
    Compute the precession matrix to transform coordinates from epoch1 to epoch2.
    """
    zeta_A, theta_A, z_A = compute_precession_angles(epoch1, epoch2)

    # Precession matrix
    P = np.dot(rotation_matrix(-z_A, 'z'), np.dot(rotation_matrix(theta_A, 'y'), rotation_matrix(-zeta_A, 'z')))

    return P


def precess_uvw(uu, v, w, ra, dec, epoch1, epoch2):
    """
    Precess the UVW coordinates from epoch1 to epoch2.
    """
    uvw = np.array([uu, v, w])

    # Compute precession matrix
    precession_matrix = compute_precession_matrix(epoch1, epoch2)

    # Apply precession matrix to UVW coordinates
    uvw_precessed = np.dot(precession_matrix, uvw)

    # Create SkyCoord object for the initial position with units
    skycoord_initial = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=FK5(equinox=epoch1))

    # Precess to the target epoch
    skycoord_target = skycoord_initial.transform_to(FK5(equinox=epoch2))

    # Calculate the phase rotation due to the coordinate change
    ra_new = skycoord_target.ra.deg
    dec_new = skycoord_target.dec.deg

    # Calculate the phase shift vector
    ra_diff = np.deg2rad(ra_new - ra)
    dec_diff = np.deg2rad(dec_new - dec)

    # Correct the UVW coordinates for the phase shift
    # Apply phase shift in RA
    uvw_corrected_ra = np.dot(rotation_matrix(ra_diff, 'z'), uvw_precessed)
    # Apply phase shift in Dec
    uvw_corrected = np.dot(rotation_matrix(dec_diff, 'x'), uvw_corrected_ra)

    return uvw_corrected, ra_new, dec_new


# Initial UVW coordinates (example values in meters)
u_initial, v_initial, w_initial = 516117.47778081, 99503.4107472, 276978.779884

# Initial epoch and target epoch
initial_epoch = Time('2021-05-13', scale='utc')
target_epoch = Time('2022-05-13', scale='utc')

# Right Ascension and Declination of the observed point (example values in degrees)
ra_initial = 242.75  # Example RA in degrees
dec_initial = 55.0  # Example Dec in degrees

# Compute precessed UVW coordinates and new RA/Dec
uvw_new, ra_new, dec_new = precess_uvw(u_initial, v_initial, w_initial, ra_initial, dec_initial, initial_epoch,
                                       target_epoch)

print(f"New UVW coordinates: u={uvw_new[0]:.3f} m, v={uvw_new[1]:.3f} m, w={uvw_new[2]:.3f} m")
print(f"Delta UVW coordinates: delta u={u_initial-uvw_new[0]:.3f} m, delta v={v_initial-uvw_new[1]:.3f} m, delta w={w_initial-uvw_new[2]:.3f} m")

print(f"New RA: {ra_new:.6f} degrees, New Dec: {dec_new:.6f} degrees")
print(f"Precession: {np.sqrt((u_initial - uvw_new[0]) ** 2 + (v_initial - uvw_new[1]) ** 2 + (w_initial - uvw_new[2]) ** 2)} meters")
