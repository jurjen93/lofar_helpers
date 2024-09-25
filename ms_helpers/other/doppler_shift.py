from casacore.tables import table
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from scipy.constants import c

def calculate_doppler_shifts(ms_path, instrument="LOFAR"):
    """
    Calculate Doppler shifts for an MS.

    Parameters:
    - ms_path (str): Path to the MS file.
    - instrument (str): Instrument name.

    Returns:
    - result: Doppler shift in kHz
    """

    # Get observation time
    with table(ms_path+"::POINTING", ack=False) as ms:
        obs_time_seconds = ms.getcol('TIME_ORIGIN')[0]  # Extract observation times in seconds
        ra, dec = ms.getcol("DIRECTION")[0][0]

    # Convert seconds to MJD
    obs_time_mjd = obs_time_seconds / 86400.0
    time = Time(obs_time_mjd, format='mjd', scale='utc')

    # Get ref frequency
    with table(ms_path+"::SPECTRAL_WINDOW", ack=False) as ms:
        observing_frequency = ms.getcol("REF_FREQUENCY")[0] * u.Hz

    # Define the location of the observatory
    observatory_location = EarthLocation.of_site(instrument)

    # Define the source coordinates
    source = SkyCoord(ra=ra * u.rad, dec=dec * u.rad, frame='icrs')

    # Calculate the radial velocity due to Earth's motion
    rv = source.radial_velocity_correction(obstime=time, location=observatory_location).to(u.m / u.s)

    # Calculate the Doppler shift for the observation
    doppler_shift = observing_frequency * (rv / c)

    return doppler_shift.to(u.kHz)
