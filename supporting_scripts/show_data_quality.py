import casacore.tables as ct
import matplotlib.pyplot as plt

t = ct.table('Abell399-401_box_1.dysco.sub.shift.avg.weights.ms.archive0')

DATA = t.getcol("DATA")
TIME = t.getcol("TIME")

#TODO: Remove?