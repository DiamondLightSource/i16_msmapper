# i16_msmapper
Simple GUI to run the msmapper code on Beamline I16

The Miller Space Mapper (msmapper) program converts x-ray diffraction scans with area detectors into reciprocal space units.
msmapper is developed by Peter Chang & SciSoft Group, Diamond Light Source Ltd.
Links:
 - [GitHub](https://github.com/DawnScience/scisoft-core/blob/master/uk.ac.diamond.scisoft.analysis/src/uk/ac/diamond/scisoft/analysis/diffraction/MillerSpaceMapper.java)
 - [Javadocs](https://alfred.diamond.ac.uk/documentation/javadocs/GDA/master/uk/ac/diamond/scisoft/analysis/diffraction/MillerSpaceMapper.html)
 - [I16 Confluence Page](https://confluence.diamond.ac.uk/display/I16/HKL+Mapping)

By Dan Porter, Diamond Light Source Ltd. 2024-2025

### Usage
```commandline
$ module load i16_msmapper
$ i16_msmapper
```

### Installation
**requirements:** *tkinter, numpy, matplotlib, h5py, hdfmap, msmapper*

**available from: https://github.com/DiamondLightSource/i16_msmapper**

Latest version from github:
```commandline
python -m pip install --upgrade git+https://github.com/DiamondLightSource/i16_msmapper.git
```

### Screenshot
![msmapper_gui](https://github.com/DiamondLightSource/i16_msmapper/blob/master/msmapper_gui.png?raw=true)


### Usage on Linux (Diamond Workstation)
The following commands can be used on a beamline or DLS linux workstation (including NXuser)
```bash
$ module load msmapper/1.9
$ python -m pip install --upgrade git+https://github.com/DiamondLightSource/i16_msmapper.git
$ python -m i16_msmapper
```

### Usage on Windows
MSMapper can be run outside Diamond by downloading the executable file. 
The following options are for Windows but files for other operating systems are available and the process is similar.

1. Install i16_msmapper as above
2. Access the MSMapper files for different operating systems here: https://alfred.diamond.ac.uk/MSMapper/master/builds-snapshot/ 
3. Download "MSMapper-1.7.0.v20240513-1606-win32.x86_64.zip" or equivalent
4. Unzip the file to your choosen location
5. Open the unzipped folder and copy the path of the executable (*shift-right-click* on `msmapperrc.exe` and select **copy as path**)
6. Run i16_msmapper using `python -m i16_msmapper`
7. From the main screen, select *Tools>set shell command*, replace the executable in the command with the copied one (**ctrl+v**)
8. Click OK
9. Use i16_msmapper as normal

# MSMapper Description
MSMapper (Miller-Space-Mapper) is a java application that re-maps detector image pixels into a 3D voxel-grid with 
coordinates defined by the reciprocal lattice and motor coordinates stored in the NeXus file. 

### MSMapper Usage
MSMapper is available on Diamond workstations via the module load system
```bash
$ module load msmapper
$ msmapper -bean /location/of/bean.json
# -OR-
$ rs_map -s 0.002 --monitor rc -o /dls/i16/data/2022/mm12345-1/processing/12345_remap.nxs /dls/i16/data/2022/mm12345-1/12345.nxs
```

### MSMapper Docs
https://alfred.diamond.ac.uk/documentation/javadocs/GDA/master/uk/ac/diamond/scisoft/analysis/diffraction/MillerSpaceMapper.html

### MSMapper Code
https://github.com/DawnScience/scisoft-core/blob/master/uk.ac.diamond.scisoft.analysis/src/uk/ac/diamond/scisoft/analysis/diffraction/MillerSpaceMapper.java


### Bean File Definition
The bean describes the options given to the MSMapper software

```text
bean = {
    inputs
      [list] List of filenames of .nxs scan files

    output
      [str] Output filename - must be in processing directory, or somewhere you can write to

    outputMode
      [str] Type of output generated, see below for options
    
    toCrystalFrame
      [bool] if False, Q-mappings will be in the basis of the lab frame, rather than the Crystal frame

    splitterName
      [str] one of the following strings "nearest", "gaussian", "negexp", "inverse"

    splitterParameter
      [float] splitter's parameter is distance to half-height of the weight function.
      If you use None or "" then it is treated as "nearest"

    scaleFactor
      [float] the oversampling factor for each image; to ensure that are no gaps in between pixels in mapping

    monitorName
      [str] sets the monitor to normalise by "rc", "ic1monitor"
    
    thirdAxis
      [list] direction of volume's third axis

    aziPlaneNormal
      [list] normal of azimuthal reference plane

    correctPolarization
      [bool] Set true to correct for polarization factor caused by transformation from laboratory frame to scattering plane

    region
      [sx, ex, sy, ey] Set rectangular region that defines the area with its bounds that contribute to output

    step
      [float, list] a single value or list if 3 values and determines the lengths of each side of the voxels in the volume

    start
      [list] location in HKL space of the bottom corner of the array.

    shape
      [list] size of the array to create for reciprocal space volume

    reduceToNonZero
      [bool] if True, attempts to reduce the volume output
}
```

### Output Modes

Q-coordinates are in the basis of the sample frame (z along c, as Busing & Levy B matrix), as defined in the input file, 
unless toCrystalFrame=False, in which case the coordinates are given in the 
lab frame (version 1.9+).

| Output Mode | Description                                                      |
|-------------|------------------------------------------------------------------|
| Area_HK     | Area in Miller-space (H,K)                                       |
| Area_KL     | Area in Miller-space (K,L)                                       |
| Area_LH     | Area in Miller-space (L,H)                                       |
| Area_QPP    | Area in q-space (parallel, perpendicular to sample surface)      |
| Area_QXY    | Area in q-space (X,Y)                                            |
| Area_QYZ    | Area in q-space (Y,Z)                                            |
| Area_QZX    | Area in q-space (Z,X)                                            |
| Coords_HKL  | Coordinates in Miller space                                      |
| Coords_Q    | Coordinates in q-space (momentum transfer)                       |
| Line_2Theta | Line in q-space (2 x theta is scattering angle, also in degrees) |
| Line_H      | Line in Miller space (H)                                         |
| Line_K      | Line in Miller space (K)                                         |
| Line_L      | Line in Miller space (L)                                         |
| Line_QX     | Line in q-space (X)                                              |
| Line_QY     | Line in q-space (Y)                                              |
| Line_QZ     | Line in q-space (Z)                                              |
| Line_Theta  | Line in q-space (2 x theta is scattering angle)                  |
| Volume_HKL  | Volume in Miller space                                           |
| Volume_Q    | Volume in q-space (crystal frame)                                |
| Volume_QCP  | Volume in cylindrical polar q-space (crystal frame)              |
| Volume_QES  | Volume in equatorial stereographic q-space (crystal frame)       |

### Bean File Example
```python
bean = {
    "inputs": ['i16/2024/mm12345-1/123456.nxs', 'i16/2024/mm12345-1/123457.nxs'],  # Filename of scan file
    "output": 'i16/2024/mm12345-1/processing/123456_rsmap.nxs',
    # Output filename - must be in processing directory, or somewhere you can write to
    "splitterName": "gaussian",  # one of the following strings "nearest", "gaussian", "negexp", "inverse"
    "splitterParameter": 2.0,
    # splitter's parameter is distance to half-height of the weight function.
    # If you use None or "" then it is treated as "nearest"
    "scaleFactor": 2.0,
    # the oversampling factor for each image; to ensure that are no gaps in between pixels in mapping
    "step": [0.002, 0.002, 0.002],
    # a single value or list if 3 values and determines the lengths of each side of the voxels in the volume
    "start": None,  # location in HKL space of the bottom corner of the array.
    "shape": None,  # size of the array to create for reciprocal space volume
    "reduceToNonZero": False  # True/False, if True, attempts to reduce the volume output
}
```