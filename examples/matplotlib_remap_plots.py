"""
matplotlib plotting or remap data

23 Sep 2024
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

file = 'data/processing/1041304_rsmap.nxs'
filename = os.path.basename(file)

with h5py.File(file, 'r') as hdf:
    # the following are links to the original scan file
    scan_command = hdf['/entry0/scan_command'].asstr()[()]  # str
    crystal = hdf['/entry0/before_scan/xtlinfo_extra/crystal_name'][()]  # str
    temp = hdf['/entry0/before_scan/temperature_controller/Tsample'][()]  # float
    unit_cell = np.reshape(hdf['/entry0/sample/unit_cell'], -1)  # 1D array, length 6
    energy = hdf['/entry0/sample/beam/incident_energy'][0]  # # float
    ubmatrix = hdf['/entry0/sample/ub_matrix'][0]  # 3D array, shape (3,3)
    # this is the processed data
    haxis = hdf['/processed/reciprocal_space/h-axis'][()]  # 1D array, length n
    kaxis = hdf['/processed/reciprocal_space/k-axis'][()]  # 1D array, length m
    laxis = hdf['/processed/reciprocal_space/l-axis'][()]  # 1D array, length o
    volume = hdf['/processed/reciprocal_space/volume'][()]  # 3D array, shape [n,m,o]

print(f"Loaded file: {filename} with volume shape: {volume.shape}")

# average angle subtended by each pixel
pixel_size = 0.172 # mm
detector_distance = 565 # mm
solid_angle = pixel_size ** 2 / detector_distance ** 2  # sr
print(f'Each pixel is normalised by the solid angle: {solid_angle: .4g} sr')

volume = volume * solid_angle

# Plot summed cuts
plt.figure(figsize=[18, 8], dpi=60)
title = f"{filename}\n{crystal} {temp:.3g} K: {scan_command}"
plt.suptitle(title, fontsize=18)

plt.subplot(131)
plt.plot(haxis, volume.sum(axis=1).sum(axis=1))
plt.xlabel('h-axis (r.l.u.)', fontsize=16)
plt.ylabel('sum axes [1,2]', fontsize=16)

plt.subplot(132)
plt.plot(kaxis, volume.sum(axis=0).sum(axis=1))
plt.xlabel('k-axis (r.l.u.)', fontsize=16)
plt.ylabel('sum axes [0,2]', fontsize=16)

plt.subplot(133)
plt.plot(laxis, volume.sum(axis=0).sum(axis=0))
plt.xlabel('l-axis (r.l.u.)', fontsize=16)
plt.ylabel('sum axes [0,1]', fontsize=16)

# Plot summed images
plt.figure(figsize=[18, 8], dpi=60)
title = f"{filename}\n{crystal} {temp:.3g} K: {scan_command}"
plt.suptitle(title, fontsize=20)
plt.subplots_adjust(wspace=0.3)

plt.subplot(131)
K, H = np.meshgrid(kaxis, haxis)
plt.pcolormesh(H, K, volume.sum(axis=2), shading='auto')
plt.xlabel('h-axis (r.l.u.)', fontsize=16)
plt.ylabel('k-axis (r.l.u.)', fontsize=16)
plt.axis('image')
#plt.colorbar()

plt.subplot(132)
L, H = np.meshgrid(laxis, haxis)
plt.pcolormesh(H, L, volume.sum(axis=1), shading='auto')
plt.xlabel('h-axis (r.l.u.)', fontsize=16)
plt.ylabel('l-axis (r.l.u.)', fontsize=16)
plt.axis('image')
#plt.colorbar()

plt.subplot(133)
L, K = np.meshgrid(laxis, kaxis)
plt.pcolormesh(K, L, volume.sum(axis=0), shading='auto')
plt.xlabel('k-axis (r.l.u.)', fontsize=16)
plt.ylabel('l-axis (r.l.u.)', fontsize=16)
plt.axis('image')
plt.colorbar()

plt.show()