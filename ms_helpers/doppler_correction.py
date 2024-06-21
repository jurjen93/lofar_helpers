import numpy as np
import casacore.tables as tables
import casacore.measures as measures
import casacore.quanta as quanta

def calculate_radial_velocity(observation_time, ant_position, source_coords, frame='LSRK'):
    me = measures.measures()
    qa = quanta.quanta()

    # Convert observation time to casacore quantity
    time_quantity = qa.quantity(observation_time, 's')

    # Set up frames
    me.doframe(me.epoch('UTC', time_quantity))
    me.doframe(me.position('ITRF', *[qa.quantity(pos, 'm') for pos in ant_position]))
    me.doframe(me.direction('J2000', *[qa.quantity(coord, 'rad') for coord in source_coords]))

    # Calculate radial velocity in the specified frame
    doppler = me.doppler('RADIO', frame)
    radial_velocity = doppler['m0']['value']  # in m/s
    return radial_velocity

def correct_doppler_shift_casacore(msname, restfreq, outms, frame='LSRK'):
    """
    Corrects Doppler shifts in a Measurement Set using casacore.

    Parameters:
    msname (str): Name of the input Measurement Set.
    restfreq (str): Rest frequency of the observed line (e.g., '1420.40575177MHz').
    outms (str): Name of the output Measurement Set with Doppler correction applied.
    frame (str): Reference frame for the Doppler correction (default is 'LSRK').

    Returns:
    None
    """
    # Initialize casacore tools
    me = measures.measures()
    qa = quanta.quanta()

    # Open the Measurement Set
    ms = tables.table(msname, readonly=False)

    # Get the observation time and antenna positions
    obs_table = tables.table(msname + '/OBSERVATION')
    obs_time = obs_table.getcol('TIME_RANGE')[0][0]
    obs_table.close()

    ant_table = tables.table(msname + '/ANTENNA')
    ant_positions = ant_table.getcol('POSITION')
    ant_table.close()

    # Get the field direction
    field_table = tables.table(msname + '/FIELD')
    phase_dir = field_table.getcol('PHASE_DIR')[0]
    field_table.close()

    # Convert the rest frequency to a casacore quantity
    restfreq_qa = qa.quantity(restfreq)

    # Prepare to calculate Doppler shifts
    me.doframe(me.epoch('UTC', qa.quantity(obs_time, 's')))
    avg_doppler = []

    for i in range(len(ant_positions)):
        ant_pos = me.position('ITRF',
                              qa.quantity(ant_positions[i][0], 'm'),
                              qa.quantity(ant_positions[i][1], 'm'),
                              qa.quantity(ant_positions[i][2], 'm'))
        me.doframe(ant_pos)
        me.doframe(me.direction('J2000',
                                qa.quantity(phase_dir[0], 'rad'),
                                qa.quantity(phase_dir[1], 'rad')))
        doppler = me.measure(me.doppler('RADIO', frame, restfreq_qa), frame)
        avg_doppler.append(doppler['m0']['value'])

    avg_doppler_shift = np.mean(avg_doppler)

    # Apply the Doppler correction factor
    doppler_factor = 1 + avg_doppler_shift / 299792.458  # speed of light in km/s

    # Update the frequencies in the SPW table
    spw_table = tables.table(msname + '/SPECTRAL_WINDOW', readonly=False)
    chan_freqs = spw_table.getcol('CHAN_FREQ')

    corrected_chan_freqs = chan_freqs / doppler_factor
    spw_table.putcol('CHAN_FREQ', corrected_chan_freqs)
    spw_table.close()

    # Copy the Measurement Set to the output file
    ms.copy(outms)
    ms.close()

    print(f"Doppler correction applied. Output MS saved as {outms}")

# Example usage:
msname = 'your_measurement_set.ms'
restfreq = '1420.40575177MHz'
outms = 'corrected_measurement_set.ms'
correct_doppler_shift_casacore(msname, restfreq, outms)