"""
This code is based on MockPyStep.py from https://dp3.readthedocs.io/en/latest/steps/PythonStep.html
Polarization conversion implemented by Jurjen de Jong.

To use this script in DPPP/DP3, you can save this file in the same folder as from where you run your DPPP command.
Or you save this python file somewhere you like and run:
export PYTHONPATH=/somewhere/you/like:$PYTHONPATH

"""

try:
    from dppp import DPStep as Step
    DP3name = 'DPPP' # default
except:
    from dp3 import Step
    DP3name = 'DP3'

from subprocess import check_output
import re
import numpy as np
import sys

#hacky way to figure out the DPPP/DP3 version (important to run this script properly)
try:
    rgx = '[0-9]+(\.[0-9]+)+'
    grep_version_string = str(check_output(DP3name+' --version', shell=True), 'utf-8')
    DP3_VERSION = float(re.search(rgx, grep_version_string).group()[0:3])
except AttributeError:
    print('WARNING: grep of DP3 version failed.')
    DP3_VERSION=0.0 # set to default

if DP3_VERSION > 5.3:
    from dp3 import Fields

class PolDiff(Step):
    """ Takes the difference of two given polarisations.
    """

    def __init__(self, parset, prefix):
        """
        Set up the step (constructor). Read the parset here.
        Set fetch_weights to True if the weights need to be read.
        Similarly for fetch_uvw.

        Args:
          parset: Parameter set for the entire pipeline
          prefix: Prefix for this step."
        """

        super().__init__()
        # try:
        #     self.pol1_str = parset.getString(prefix + "pol1")
        # except RuntimeError:
        #     self.pol1_str = False
        # except AttributeError:
        #     # DP3 Python bindings have been renamed.
        #     try:
        #         self.pol1_str = parset.get_string(prefix + "pol1")
        #     except RuntimeError:
        #         self.pol1_str = False

        # try:
        #     self.pol2_str = parset.get_string(prefix + "pol2")
        # except RuntimeError:
        #     self.pol2_str=False
        # except AttributeError:
        #     # DP3 Python bindings have been renamed.
        #     try:
        #         self.pol2_str = parset.get_string(prefix + "pol2")
        #     except RuntimeError:
        #         self.pol2_str=False

        # pols_lin = {'XX':0, 'XY':1, 'YX':2, 'YY':3}
        # pols_circ = {'RR':0, 'RL':1, 'LR':2, 'LL':3}

        # if self.pol1_str in pols_lin.keys():
        #     self.pol1 = pols_lin[self.pol1_str]
        # elif self.pol1_str in pols_circ.keys():
        #     self.pol1 = pols_circ[self.pol1_str]
        # else:
        #     sys.exit('Unknown polarisation specified for pol1.')

        # if self.pol2_str in pols_lin.keys():
        #     self.pol2 = pols_lin[self.pol2_str]
        # elif self.pol2_str in pols_circ.keys():
        #     self.pol2 = pols_circ[self.pol2_str]
        # else:
        #     sys.exit('Unknown polarisation specified for pol2.')
        self.fetch_uvw = True

    def get_required_fields(self):
        if DP3_VERSION>5.3:
            return (Fields.DATA | Fields.FLAGS | Fields.WEIGHTS | Fields.UVW)
        else:
            pass

    def get_provided_fields(self):
        if DP3_VERSION>5.3:
            return Fields()
        else:
            pass

    def update_info(self, dpinfo):
        """
        Process metadata. This will be called before any call to process.

        Args:
          dpinfo: DPInfo object with all metadata, see docs in pydp3.cc
        """

        super().update_info(dpinfo)

        # Make sure data is read
        self.info().set_need_vis_data()

        # Make sure data will be written
        self.info().set_write_data()

    def show(self):
        """Print a summary of the step and its settings"""

        print("\nPolDiff")
        print("\nCreating XX-YY phase difference.")
        # print("\nSubtracting {:s} from {:s}".format(self.pol1_str, self.pol2_str))


    def process(self, dpbuffer):
        """
        Process one time slot of data. This function MUST call process_next_step.

        Args:
          dpbuffer: DPBuffer object which can contain data, flags and weights
                    for one time slot.
        """

        data = np.array(dpbuffer.get_data(), copy=False, dtype=np.complex64)
        # newdata = data[:, :, self.pol1] - data[:, :, self.pol2]
        phasediff =  np.copy(np.angle(data[:, :, 0]) - np.angle(data[:, :, 3]))  #RR - LL
        data[:, :, 0] = 0.5 * np.exp(1j * phasediff)  # because I = RR+LL/2 (this is tricky because we work with phase diff)
        data[:, :, 3] = data[:, :, 0]

        # Send processed data to the next step
        if DP3_VERSION>5.3:
            next_step = self.get_next_step()
            if next_step is not None:
                next_step.process(dpbuffer)
        else:
            self.process_next_step(dpbuffer)

    def finish(self):
        """
        If there is any remaining data, process it. This can be useful if the
        step accumulates multiple time slots.
        """
        print("\nCreated XX-YY phase difference.")
        # print("\nSubtracted {:s} from {:s}".format(self.pol1_str, self.pol2_str))