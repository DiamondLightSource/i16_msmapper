"""
matplotlib voxel plot for 3D plotting of remapped data

3D voxel / volumetric plot
https://matplotlib.org/stable/gallery/mplot3d/voxels.html
https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.voxels.html

23 Sep 2024
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import time

file = 'data/processing/1041304_rsmap.nxs'
file = '/dls/i16/data/2025/cm40634-12/processing/1105476_workflows_msmapper.nxs'
filename = os.path.basename(file)

with h5py.File(file, 'r') as hdf:
    # the following are links to the original scan file
    scan_command = hdf['/entry0/scan_command'].asstr()[()]  # str
    # crystal = hdf['/entry0/before_scan/xtlinfo_extra/crystal_name'][()]  # str
    # temp = hdf['/entry0/before_scan/temperature_controller/Tsample'][()]  # float
    crystal = hdf['/entry0/sample/name'][()]  # str
    temp = hdf['/entry0/instrument/temperature_controller/Tsample'][()]  # float
    unit_cell = np.reshape(hdf['/entry0/sample/unit_cell'], -1)  # 1D array, length 6
    energy = hdf['/entry0/sample/beam/incident_energy'][...]  # # float
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

# rebin to make smaller plot
# haxis = haxis[::2]
# kaxis = kaxis[::2]
# laxis = laxis[::2]
# volume = volume[::2, ::2, ::2]

# Coordinates - plt.voxel takes coordinates of each corner, so hkl must be extended or volume reduced
kk, hh, ll = np.meshgrid(kaxis, haxis, laxis)
# volume = volume[:-1, :-1, :-1]
print(f"h-axis: {hh.shape}, min={haxis.min()}, max={haxis.max()}")
print(f"k-axis: {kk.shape}, min={kaxis.min()}, max={kaxis.max()}")
print(f"l-axis: {ll.shape}, min={laxis.min()}, max={laxis.max()}")
print(f"volume: {volume.shape}, min={volume.min()}, max={volume.max()}")


title = f"{filename}\n{crystal} {temp:.3g} K: {scan_command}"
cmap = plt.get_cmap('Greys')
min_intensities = [0.001, 0.01, 0.1, 0.2]
# alphas = np.arange(0.1, 1, 0.9 / len(min_intensities))
alphas = np.linspace(0.1, 1, len(min_intensities))
max_volume = volume.max()

# Plot 3D point cloud
fig = plt.figure(figsize=(12, 10), dpi=100)
ax = fig.add_subplot(projection='3d')
t0 = time.time()
for cut, alpha in zip(min_intensities, alphas):
    t1 = time.time()
    filled = volume > (cut * max_volume)
    nfilled = np.sum(filled)
    # print(f"cut={cut}, filled = {nfilled}({nfilled/volume.size:.3%})")
    color = cmap(alpha, alpha=alpha)
    ax.plot(hh[filled].flatten(), kk[filled].flatten(), ll[filled].flatten(), '.', c=color)
    t2 = time.time()
    # print(f"cut time = {t2 - t1:.2f}s")
t2 = time.time()
print(f"Point cloud total time: {t2 - t0:.2f}s")

# Plot 3D voxel volumetrix plot
fig = plt.figure(figsize=(12, 10), dpi=100)
ax = fig.add_subplot(projection='3d')
t0 = time.time()
for cut, alpha in zip(min_intensities, alphas):
    t1 = time.time()
    filled = volume > (cut * max_volume)
    nfilled = np.sum(filled)
    color = cmap(alpha, alpha=alpha)
    # rebin transparent voxels for speed
    if cut < 0.1:
        ax.voxels(
            hh[::2, ::2, ::2],
            kk[::2, ::2, ::2],
            ll[::2, ::2, ::2],
            filled=filled[::2, ::2, ::2][:-1, :-1, :-1],
            facecolors=color
        )
    else:
        ax.voxels(hh, kk, ll, filled[:-1, :-1, :-1], facecolors=color)
    t2 = time.time()
    print(f"voxel cut={cut}, filled = {nfilled}({nfilled/volume.size:.3%}), time: {t2 - t1:.2f}s")
t2 = time.time()
print(f"Voxel total time: {t2 - t0:.2f}s")

ax.set_xlabel('h-axis')
ax.set_ylabel('k-axis')
ax.set_zlabel('l-axis')
ax.set_title(title)

# Plot histogram
ax = plt.figure().add_subplot()
n, bins, patches = ax.hist(np.log10(volume[volume > 0].flatten()), 100)

for cut, alpha in zip(min_intensities, alphas):
    logval = np.log10(cut * max_volume)
    ax.axvline(logval, c='k')
    for patch in patches:
        if patch.xy[0] >= logval:
            patch.set_color(cmap(alpha))
# ax.loglog()
ax.set_xlabel('Log$_{10}$ Voxel Intensity')
ax.set_ylabel('N')
ax.set_title(title)

plt.show()
