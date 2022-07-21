# Based on MockPyStep.py from https://dp3.readthedocs.io/en/latest/steps/PythonStep.html

#export PYTHONPATH=/path/to/dp3/lib64/python3.6/site-packages:$PYTHONPATH

import parameterset
from pydp3 import Step
import numpy as np

class DataConv(Step):
    """
    Convert data.

    datafactor --> multiply data with a specific factor
    weightfactor --> multiply weights with a specific factor

    lin2circ --> convert from linear to circular data
    circ2lin --> convert from circular to linear data
    """
    def __init__(self, parset, prefix):
        """
        Set up the step (constructor). Read the parset here.
        Set fetch_weights to True if the weights need to be read.
        Similarly for fetch_uvw.

        Args:
          parset: Parameter set for the entire pipeline
          prefix: Prefix for this step, e.g. "thisstepname."
        """

        super().__init__()
        self.datafactor = parset.getDouble(prefix + "datafactor")
        self.weightsfactor = parset.getDouble(prefix + "weightsfactor")
        self.lin2circ = parset.getDouble(prefix + "lin2circ")
        self.circ2lin = parset.getDouble(prefix + "circ2lin")

        self.fetch_weights = True
        self.fetch_uvw = False

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
        print("\nMockPyStep")
        print(f"  data factor:    {self.datafactor}")
        print(f"  weights factor: {self.weightsfactor}")

    def process(self, dpbuffer):
        """
        Process one time slot of data. This function MUST call process_next_step.

        Args:
          dpbuffer: DPBuffer object which can contain data, flags and weights
                    for one time slot.
        """

        data = np.array(dpbuffer.get_data(), copy=False)
        weights = np.array(dpbuffer.get_weights(), copy=False)

        # Do the operation on data

        if self.circ2lin:

            """
            circ2lin
            XX =   RR +  RL +  LR +  LL
            XY =  iRR - iRL + iLR - iLL
            YX = -iRR - iRL + iLR + iLL
            YY =   RR -  RL -  LR +  LL
            """

            data = np.transpose(np.array([
                0.5 * (data[:, :, 0] + data[:, :, 1] + data[:, :, 2] + data[:, :, 3]),
                0.5 * (1j * data[:, :, 0] - 1j * data[:, :, 1] + 1j * data[:, :, 2] - 1j * data[:, :, 3]),
                0.5 * (-1j * data[:, :, 0] - 1j * data[:, :, 1] + 1j * data[:, :, 2] + 1j * data[:, :, 3]),
                0.5 * (data[:, :, 0] - data[:, :, 1] - data[:, :, 2] + data[:, :, 3])]),
                (1, 2, 0))

        elif self.lin2circ:
            """
            lin2circ
            RR = XX - iXY + iYX + YY
            RL = XX + iXY + iYX - YY
            LR = XX - iXY - iYX - YY
            LL = XX + iXY - iYX + YY
            """

            data = np.transpose(np.array([
                0.5 * (data[:, :, 0] - 1j * data[:, :, 1] + 1j * data[:, :, 2] + data[:, :, 3]),
                0.5 * (data[:, :, 0] + 1j * data[:, :, 1] + 1j * data[:, :, 2] - data[:, :, 3]),
                0.5 * (data[:, :, 0] - 1j * data[:, :, 1] - 1j * data[:, :, 2] - data[:, :, 3]),
                0.5 * (data[:, :, 0] + 1j * data[:, :, 1] - 1j * data[:, :, 2] + data[:, :, 3])]),
                (1, 2, 0))

        data *= self.datafactor
        weights *= self.weightsfactor

        # Send processed data to the next step
        self.process_next_step(dpbuffer)

    def finish(self):
        """
        If there is any remaining data, process it. This can be useful if the
        step accumulates multiple time slots.
        """
        pass