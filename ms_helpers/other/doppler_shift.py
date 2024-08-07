from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from scipy.constants import c

# Define the observing frequency
observing_frequency = 140 * u.MHz

# Define the location of the observatory (example: Very Large Array)
observatory_location = EarthLocation.of_site('LOFAR')

# Define the source coordinates
source = SkyCoord(ra=-2.04639855*u.rad, dec=0.95905842*u.rad, frame='icrs')

for t in ['2018-11-26', '2020-05-24', '2020-11-14', '2021-5-13']:
    print(t)

    # Define the observation times
    time = Time(t)

    # Calculate the AltAz frame for the observation times
    altaz = AltAz(location=observatory_location, obstime=time)

    # Transform the source coordinates to the AltAz frame
    altaz_coord = source.transform_to(altaz)

    # Calculate the radial velocity due to Earth's motion
    rv = altaz_coord.radial_velocity_correction()

    # Calculate the Doppler shift for each observation
    doppler_shift = observing_frequency * (rv / c / u.m*u.s)

    print(f"Doppler shift for observation 1: {doppler_shift.to(u.kHz)}")

    # Adjust the observed frequencies before stacking
    adjusted_frequency = observing_frequency + doppler_shift

    print(f"Adjusted frequency for observation 1: {adjusted_frequency.to(u.MHz)}")
