#!/usr/bin/env python3
import matplotlib.pyplot as plt
from nptdms import TdmsFile
from prefixed import Float
import cavity_pert


cr_001_a = TdmsFile.read("./CR-001-a.tdms")
empty_sample = TdmsFile.read("./Empty Sample.tdms")

results = cavity_pert.analyse(cr_001_a, empty_sample, cavity_pert.CAVITY_SQUARED_RADII["Gandalf"], 0.001**2)

for name, data in results.items():
    data.plot(f"{name} = {Float(data.mean_f):.2h} Hz")

plt.show()
