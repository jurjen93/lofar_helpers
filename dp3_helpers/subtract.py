"""
!!!DOES NOT WORK (YET) !!!
"""

import numpy as np
from dp3 import Fields
from dp3.pydp3 import Step

class SubtractModelData(Step):
    """DPStep to subtract MODEL_DATA from DATA"""

    def __init__(self, parset, prefix):
        super().__init__()
        self.model_data = parset.get_string(prefix + "model_data")

    def show(self):
        print("\nDATA - MODEL_DATA")

    def get_required_fields(self):
        return Fields()

    def get_provided_fields(self):
        return Fields()

    def process(self, dpbuffer):

        # input data
        model_data = np.array(dpbuffer.get_data(self.model_data), copy=False)
        data = np.array(dpbuffer.get_data(), copy=False)

        # subtract
        data -= model_data
        self.get_next_step().process(dpbuffer)

    def finish(self):
        pass