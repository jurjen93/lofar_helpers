from __future__ import print_function, division
import numpy as np
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import (SkyCoord, FK5, CIRS, EarthLocation, CartesianRepresentation,
                                 solar_system_ephemeris, get_body_barycentric_posvel)

def calculate_uvw(X1, Y1, Z1, X2, Y2, Z2, time_str, alpha, delta):
    """
    Calculate UVW coordinates for a given observation.

    Parameters:
    X1, Y1, Z1 : float : Geocentric coordinates of antenna 1
    X2, Y2, Z2 : float : Geocentric coordinates of antenna 2
    time_str : str : Observation time in UTC (ISO format string)
    alpha : float : Right Ascension of the source in degrees
    delta : float : Declination of the source in degrees
    latitude : float : Latitude of the observation location in degrees
    longitude : float : Longitude of the observation location in degrees

    Returns:
    uvw : numpy.ndarray : UVW coordinates
    """

    def midpoint_to_lat_lon(X1, Y1, Z1, X2, Y2, Z2):
        # Calculate midpoint coordinates
        X_m = (X1 + X2) / 2
        Y_m = (Y1 + Y2) / 2
        Z_m = (Z1 + Z2) / 2

        # Create a Cartesian representation of the midpoint
        mid_cartesian = CartesianRepresentation(X_m * u.m, Y_m * u.m, Z_m * u.m)

        # Convert to EarthLocation to get latitude and longitude
        mid_location = EarthLocation.from_geocentric(mid_cartesian.x, mid_cartesian.y, mid_cartesian.z)

        return mid_location.lat.deg, mid_location.lon.deg

    # Baseline vector
    B = np.array([X2 - X1, Y2 - Y1, Z2 - Z1])

    # Define observer location
    latitude, longitude = midpoint_to_lat_lon(X1, Y1, Z1, X2, Y2, Z2)
    location = EarthLocation(lat=latitude, lon=longitude)

    # Convert to LST
    observation_time = Time(time_str, location=location)
    lst = observation_time.sidereal_time('mean', longitude=longitude).deg

    # Calculate the hour angle
    H = lst - alpha
    if H < 0:
        H += 360.0

    # Convert degrees to radians
    H_rad = np.radians(H)
    delta_rad = np.radians(delta)

    # UVW Transformation matrix
    transformation_matrix = np.array([
        [np.sin(H_rad), np.cos(H_rad), 0],
        [-np.sin(delta_rad) * np.cos(H_rad), np.sin(delta_rad) * np.sin(H_rad), np.cos(delta_rad)],
        [np.cos(delta_rad) * np.cos(H_rad), -np.cos(delta_rad) * np.sin(H_rad), np.sin(delta_rad)]
    ])

    # Compute UVW coordinates
    uvw = np.dot(transformation_matrix, B)
    return uvw


# def aberration_shift(timestamp, ra, dec):
#
#     """ FROM: https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/aberration.html """
#
#     # Convert calendar date to JD using the datetime package
#     dt = datetime.datetime.strptime(timestamp, "%Y-%m-%dT%H:%M:%S")
#     jd = pyasl.jdcnv(dt)
#
#     # Specify RA and DEC
#     ra = np.atleast_1d(ra)
#     dec = np.atleast_1d(dec)
#
#     # Get change in RA and DEC due to annual aberration
#     res = pyasl.co_aberration(np.repeat(jd, ra.size), ra, dec)
#
#     return res[0][0]*3600, res[1][0]*3600


def aberration_correction(timestamp, ra, dec):
    """
    Calculate the aberration corrected RA and DEC for given input RA, DEC and timestamp.

    Parameters:
    ra (float): Right Ascension in degrees
    dec (float): Declination in degrees
    timestamp (str): Timestamp in ISO format (e.g., '2024-07-17T00:00:00')

    Returns:
    corrected_ra (float): Aberration corrected Right Ascension in degrees
    corrected_dec (float): Aberration corrected Declination in degrees
    """
    # Convert timestamp to Julian Date
    time = Time(timestamp, format='isot', scale='utc')

    # Initial coordinates
    # sky_coord = SkyCoord(ra=ra, dec=dec, unit='deg')

    # Calculate the position and velocity of Earth
    with solar_system_ephemeris.set('builtin'):
        pos, vel = get_body_barycentric_posvel('earth', time)

    # Speed of light in AU/day
    c = 173.144632674240

    # Aberration correction
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    # Earth's velocity components
    vx = vel.x.value
    vy = vel.y.value
    vz = vel.z.value

    # Proper motion in RA and DEC
    delta_ra = (-vx * np.sin(ra_rad) + vy * np.cos(ra_rad)) / c
    delta_dec = (-vx * np.cos(ra_rad) * np.sin(dec_rad)
                 - vy * np.sin(ra_rad) * np.sin(dec_rad)
                 + vz * np.cos(dec_rad)) / c


    return np.rad2deg(delta_ra) * 3600, np.rad2deg(delta_dec) * 3600


def precession_nutation_shift(timestamp, ra, dec):
    """
    Calculate the precession and nutation corrected RA and DEC for given input RA, DEC, and timestamp.

    Parameters:
    ra (float): Right Ascension in degrees
    dec (float): Declination in degrees
    timestamp (str): Timestamp in ISO format (e.g., '2024-07-17T00:00:00')

    Returns:
    corrected_ra (float): Precession and nutation corrected Right Ascension in degrees
    corrected_dec (float): Precession and nutation corrected Declination in degrees
    """
    # Convert timestamp to astropy Time object
    t = Time(timestamp)

    # Create a SkyCoord object for the input RA and DEC
    coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')

    # Convert the coordinates to the FK5 frame at the given epoch to account for precession
    coord_fk5 = coord.transform_to(FK5(equinox=t))

    # Convert the coordinates to the CIRS frame to account for nutation
    coord_cirs = coord_fk5.transform_to(CIRS(obstime=t))

    # Extract the corrected RA and DEC
    corrected_ra = coord_cirs.ra.degree
    corrected_dec = coord_cirs.dec.degree

    return (ra-corrected_ra)*3600, (dec-corrected_dec)*3600


def calculate_uvw_shift(ra_deg, dec_deg, baseline_length, ra_shift_arcsec, dec_shift_arcsec):
    """
    Calculate the shift in UVW coordinates based on shifts in RA and DEC.

    Parameters:
    ra_deg (float): Original Right Ascension in degrees
    dec_deg (float): Original Declination in degrees
    baseline_length (float): Baseline length in meters
    ra_shift_arcsec (float): Shift in Right Ascension in arcseconds
    dec_shift_arcsec (float): Shift in Declination in arcseconds

    Returns:
    delta_u (float): Shift in U coordinate in meters
    delta_v (float): Shift in V coordinate in meters
    delta_w (float): Shift in W coordinate in meters
    """

    def ra_dec_to_radians(ra_shift_arcsec, dec_shift_arcsec):
        """Convert RA and DEC shifts from arcseconds to radians."""
        ra_shift_rad = np.deg2rad(ra_shift_arcsec / 3600.0)
        dec_shift_rad = np.deg2rad(dec_shift_arcsec / 3600.0)
        return ra_shift_rad, dec_shift_rad

    def direction_cosines(ra_rad, dec_rad):
        """Calculate direction cosines from RA and DEC in radians."""
        l = np.cos(dec_rad) * np.cos(ra_rad)
        m = np.cos(dec_rad) * np.sin(ra_rad)
        n = np.sin(dec_rad)
        return l, m, n

    def uvw_shift(ra_rad, dec_rad, ra_shift_rad, dec_shift_rad, baseline_length):
        """Calculate the shift in UVW coordinates based on RA and DEC shifts."""
        # Original direction cosines
        l0, m0, n0 = direction_cosines(ra_rad, dec_rad)

        # Shifted RA and DEC
        ra_shifted_rad = ra_rad + ra_shift_rad
        dec_shifted_rad = dec_rad + dec_shift_rad

        # Shifted direction cosines
        l1, m1, n1 = direction_cosines(ra_shifted_rad, dec_shifted_rad)

        # Difference in direction cosines
        delta_l = l1 - l0
        delta_m = m1 - m0
        delta_n = n1 - n0

        # Convert direction cosine shifts to UVW coordinate shifts
        delta_u = baseline_length * delta_l
        delta_v = baseline_length * delta_m
        delta_w = baseline_length * delta_n

        return delta_u, delta_v, delta_w

    # Convert shifts to radians
    ra_shift_rad, dec_shift_rad = ra_dec_to_radians(ra_shift_arcsec, dec_shift_arcsec)

    # Convert RA and DEC to radians
    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)

    # Calculate UVW coordinate shifts
    delta_u, delta_v, delta_w = uvw_shift(ra_rad, dec_rad, ra_shift_rad, dec_shift_rad, baseline_length)

    return delta_u, delta_v, delta_w


# Example usage:
X1, Y1, Z1 = 3826896.235000, 460979.455000, 5064658.203000
X2, Y2, Z2 = 3370271.657000, 712125.881000, 5349991.165000
ra = 241.88  # Right Ascension in degrees
dec = 54.33   # Declination in degrees

# Length of a sidereal day
sidereal_day = TimeDelta(23.9344696 * 3600, format='sec')  # 23 hours, 56 minutes, 4 seconds


for n, timez in enumerate(["2018-11-26T07:13:43", "2020-05-24T19:20:26", "2020-11-14T08:11:00", "2021-05-13T19:41:00"]):

    dra_aber, ddec_aber = aberration_correction(timez, ra, dec)
    dra_prec, ddec_prec = precession_nutation_shift(timez, ra, dec)
    # if n>=1:
    #     ddec_aber-=40
    du, dv, dw = calculate_uvw_shift(ra, dec, 600_000, dra_prec + dra_aber, ddec_prec + ddec_aber)
    print(timez)
    print(dra_prec, ddec_prec)
    print(dra_aber, ddec_aber)


    if n==0:
        ref_dra_prec, ref_ddec_prec = dra_prec, ddec_prec
        ref_dra_aber, ref_ddec_aber = dra_aber, ddec_aber
        ref_du, ref_dv, ref_dw = du, dv, dw

    # else:
    #     print(dra_prec-ref_dra_prec, ddec_prec-ref_ddec_prec)
    #     print(dra_aber-ref_dra_aber, ddec_aber-ref_ddec_aber)
    #     print((dra_prec-ref_dra_prec)-(dra_aber-ref_dra_aber), (ddec_prec-ref_ddec_prec)-(ddec_aber-ref_ddec_aber))



    print(du-ref_du, dv-ref_dv)
    if n==0:
        plt.scatter([(du-ref_du)/1000+518.45-i*0.03 for i in range(10)],
                    [(dv-ref_dv)/1000+212.55 + 0.3*i for i in range(10)], label=timez, marker="*")
    else:
        plt.scatter([(du-ref_du)/1000-0.033+518.45-i*0.03 for i in range(10)],
                    [(dv-ref_dv)/1000-0.033+212.55 + 0.3*i for i in range(10)], label=timez, marker="*")
    print()

    # # pos = np.array(
    # #     [calculate_uvw(X1, Y1, Z1, X2, Y2, Z2,
    # #                    (Time(timez) + i*4/(23.9344696 * 3600)).iso,
    # #                    ra-(dra_prec+dra_aber)/3600,
    # #                    dec-(ddec_prec+ddec_aber)/3600)
    # #      for i in range(int(3600))])
    # pos = np.array(
    #     [calculate_uvw(X1, Y1, Z1, X2, Y2, Z2,
    #                    (Time(timez) + i*8/(23.9344696 * 3600)).iso,
    #                    ra+(dra_prec+dra_aber)/3600,
    #                    dec+(ddec_prec+ddec_aber)/3600)
    #      for i in range(int(3600))])
    # plt.scatter(pos[:, 0]/1000, pos[:, 1]/1000, label=timez, s=30, marker="*")

plt.xlim(518.1, 518.6)
plt.ylim(212.4, 215.1)
plt.legend()
plt.show()
