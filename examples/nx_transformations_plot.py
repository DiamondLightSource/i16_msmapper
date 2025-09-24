"""
Example Script

September 2025
"""

import matplotlib.pyplot as plt
import h5py
from i16_msmapper.nx_transformations import NXScan

with h5py.File('data/1109527.nxs') as nxs:
    scan = NXScan(nxs)
    scan.plot_wavevectors()

    scan.plot_instrument()

    plt.show()

