"""
Runs msmapper on Diamond linux workstations

run with
$ module load msmapper

"""

import sys
import os
import tempfile
import subprocess
import json
from typing import TypedDict
from typing_extensions import Unpack

import numpy as np
import hdfmap

SHELL_CMD = "msmapper -bean %s"
TEMP_BEAN = 'tmp_remap.json'
TEMP_NEXUS = 'tmp_remap.nxs'
TEMPLATE = os.path.abspath(os.path.join(os.path.dirname(__file__), 'templates', 'msmapper_script_template.txt'))
PLOT_TEMPLATE = os.path.abspath(os.path.join(os.path.dirname(__file__), 'templates', 'plotter_script_template.txt'))

# Find writable directory
TEMPDIR = tempfile.gettempdir()
if not os.access(TEMPDIR, os.W_OK):
    TEMPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if not os.access(TEMPDIR, os.W_OK):
        TEMPDIR = os.path.expanduser('~')
print('Writable TEMPDIR = %s' % TEMPDIR)
TEMP_BEAN = os.path.join(TEMPDIR, TEMP_BEAN)
TEMP_NEXUS = os.path.join(TEMPDIR, TEMP_NEXUS)


class Options(TypedDict):
    input_files: list[str]
    output_file: str
    start: None | list[float, float, float] = None
    shape: None | list[int, int, int] = None
    step: None | list[float] | list[float, float, float] = None
    output_mode: None | str = None
    normalisation: None | str = None
    polarisation: None | bool = None
    detector_region: None | list[int, int, int, int] = None
    reduce_box: None | bool = None
    third_axis: None | list[float, float, float] = None
    azi_plane_normal: None | list[float, float, float] = None


def msmapper_version():
    """Return numeric version of msmapper, or 0 if not available"""
    path = next((p for p in sys.path if 'apps/msmapper/' in p), '')  # find /dls_sw/apps/msmapper/version
    print('path: ', path)
    try:
        return float(os.path.basename(path))
    except ValueError:
        return 0


def msmapper(bean_file):
    """
    Run msmapper in subprocess, requires to be in msmapper module
      (python 3.9)$ msmapper -bean bean_file
    :param bean_file: str location of json file with input options
    :return: Returns on completion
    """
    print('\n\n\n################# Starting msmapper ###################')
    print(f"Running command:\n{SHELL_CMD % bean_file}\n\n\n")
    output = subprocess.run(SHELL_CMD % bean_file, shell=True, capture_output=True)
    print(output.stdout.decode())
    print('\n\n\n################# msmapper finished ###################\n\n\n')


def batch_commands(cmd_list):
    """
    Run a sequence of commnds in the terminal
    :param cmd_list: list of commands
    :return: Returns on completion
    """
    print('\n\n\n################# Starting msmapper ###################\n\n\n')
    commands = '\n'.join(cmd_list)
    output = subprocess.run(commands, shell=True, capture_output=True)
    print(output.stdout.decode())
    print('\n\n\n################# msmapper finished ###################\n\n\n')


def inspect_file(nexus_file):
    """
    Inspect NeXus file and return string output
    :param nexus_file: str filename of NeXus file
    :return: str
    """
    scan = hdfmap.NexusLoader(nexus_file)
    return scan.summary()


def get_nexus_bmatrix(nexus_file):
    """
    Get the B-matrix using unit cell parameters form nexus file
    Returns the B-Matrix in inverse Angstroms (units of 2pi) in the Busing & Levy formalism
    """
    scan = hdfmap.NexusLoader(nexus_file)
    if 'unit_cell' not in scan.map:
        raise KeyError('Unit Cell parameters are not in nexus file')

    a, b, c, alpha, beta, gamma = scan('unit_cell')

    alpha1 = np.deg2rad(alpha)
    alpha2 = np.deg2rad(beta)
    alpha3 = np.deg2rad(gamma)

    beta1 = np.arccos((np.cos(alpha2) * np.cos(alpha3) - np.cos(alpha1)) / (np.sin(alpha2) * np.sin(alpha3)))
    beta2 = np.arccos((np.cos(alpha1) * np.cos(alpha3) - np.cos(alpha2)) / (np.sin(alpha1) * np.sin(alpha3)))
    beta3 = np.arccos((np.cos(alpha1) * np.cos(alpha2) - np.cos(alpha3)) / (np.sin(alpha1) * np.sin(alpha2)))

    b1 = 1 / (a * np.sin(alpha2) * np.sin(beta3))
    b2 = 1 / (b * np.sin(alpha3) * np.sin(beta1))
    b3 = 1 / (c * np.sin(alpha1) * np.sin(beta2))

    # c1 = b1 * b2 * np.cos(beta3)
    # c2 = b1 * b3 * np.cos(beta2)
    # c3 = b2 * b3 * np.cos(beta1)

    b_matrix = np.array([
        [b1, b2 * np.cos(beta3), b3 * np.cos(beta2)],
        [0, b2 * np.sin(beta3), -b3 * np.sin(beta2) * np.cos(alpha1)],
        [0, 0, 1 / c]
    ])
    return 2 * np.pi * b_matrix


def get_nexus_hkl(nexus_file):
    """
    Get hkl values from Nexus file
    """
    scan = hdfmap.NexusLoader(nexus_file)
    if 'h_axis' in scan.map:
        h, k, l = scan('h_axis.mean(), k_axis.mean(), l_axis.mean()')
    elif 'diffractometer_sample_h' in scan.map:
        h, k, l = scan.get_data(*[
            'diffractometer_sample_h',
            'diffractometer_sample_k',
            'diffractometer_sample_l',
        ])
    elif 'h' in scan.map:
        print(f"Found 'h' axis at: {scan.map.get_path('h')}")
        h, k, l = scan.get_data(*['h', 'k', 'l'])
    else:
        raise KeyError('h,k,l are not in the nexus file')
    return h, k, l


def get_remap_values(nexus_file):
    """
    Get hkl values from Nexus file
    """
    scan = hdfmap.NexusLoader(nexus_file)
    if 'h_axis' in scan.map:
        h_axis, k_axis, l_axis = scan('h_axis, k_axis, l_axis')
    else:
        raise KeyError('h,k,l are not in the nexus file')
    h_min = h_axis[0]
    h_step = np.mean(np.diff(h_axis))
    h_size = len(h_axis)
    k_min = k_axis[0]
    k_step = np.mean(np.diff(k_axis))
    k_size = len(k_axis)
    l_min = l_axis[0]
    l_step = np.mean(np.diff(l_axis))
    l_size = len(l_axis)
    return h_min, h_step, h_size, k_min, k_step, k_size, l_min, l_step, l_size


def get_pixel_steps(nexus_file, nxs_file=TEMP_NEXUS, bean_file=TEMP_BEAN):
    """
    Get minimum pixel steps for a scan file
     - Creates & stores json bean file with outputMode: 'Coords_HKL' + fixed pixel indexes
     - Runs msmapper (requires msmapper environment)
     - Opens resulting Nexus file and loads coordinates
     - determines variation in hkl between adjacent pixels

    :param nexus_file: str filename of .nxs scan file
    :return: h_diff, k_diff, l_diff
    """

    bean = {
        "inputs": [nexus_file],
        "output": nxs_file,
        "outputMode": "Coords_HKL",
        "pixelIndexes": [[0, 0, 0], [1, 1, 1], [2, 2, 2]],
    }
    json.dump(bean, open(bean_file, 'w'), indent=4)
    print('bean file written to: %s' % bean_file)

    msmapper(bean_file)

    # 3. Read nxs file
    print('\nReading %s' % nxs_file)
    coords = hdfmap.hdf_data(nxs_file, '/processed/reciprocal_space/coordinates')

    hkl_diff = np.abs(np.mean(np.diff(coords, axis=0), axis=0))
    print('\n\n***Results***')
    print('File: %s' % nexus_file)
    print('pixel step (h,k,l) = (%.2g, %.2g, %.2g)\n\n' % (hkl_diff[0], hkl_diff[1], hkl_diff[2]))
    return hkl_diff


def generate_pixel_coordinates(nexus_file, bean_file=TEMP_BEAN):
    """
    Get minimum pixel steps for a scan file
     - Creates & stores json bean file with outputMode: 'Coords_HKL' + fixed pixel indexes
     - Runs msmapper (requires msmapper environment)
     - Opens resulting Nexus file and loads coordinates
     - determines variation in hkl between adjacent pixels

    :param nexus_file: str filename of .nxs scan file
    :return: nexus_filename
    """

    nxs_map = hdfmap.create_nexus_map(nexus_file)
    path = nxs_map.get_image_path()
    if not path:
        raise Exception(f"No detector found in file: {nexus_file}")
    det_shape = nxs_map.datasets[path].shape

    pixelIndexes = []
    for idx in range(det_shape[0]):
        pixelIndexes.append([idx, 0, 0])
        pixelIndexes.append([idx, det_shape[1], 0])
        pixelIndexes.append([idx, det_shape[1], det_shape[2]])
        pixelIndexes.append([idx, 0, det_shape[2]])
        pixelIndexes.append([idx, 0, 0])

    outfile = nexus_file.replace('.nxs', '_pixel_hkl.nxs')
    bean = {
        "inputs": [nexus_file],
        "output": outfile,
        "outputMode": "Coords_HKL",
        "pixelIndexes": pixelIndexes,
    }
    json.dump(bean, open(bean_file, 'w'), indent=4)
    print('bean file written to: %s' % bean_file)

    msmapper(bean_file)

    # 3. Read nxs file
    print('\nReading %s' % outfile)
    coords = hdfmap.hdf_data(outfile, '/processed/reciprocal_space/coordinates')
    return coords


def rsmap_command(input_files, output_file, step=0.002):
    """
    Create rsmap command from values
    $ rs_map -s 0.002 -o /dls/i16/data/2022/mm12345-1/processing/12345_remap.nxs /dls/i16/data/2022/mm12345-1/12345.nxs
    :param input_files: list of scan file locations
    :param output_file: str localtion of output file
    :param step: [dh, dk, dl] step size in each direction - size of voxel in reciprocal lattice units
    :return: str command
    """
    input_files = np.asarray(input_files, dtype=str).reshape(-1).tolist()
    return f"rs_map -s {step} -o {output_file} {input_files[0]}"


def rsmap_batch(input_files, output_directory, step=0.002):
    """
    Create batch of rsmap command for many files
    $ rs_map -s 0.002 -o /dls/i16/data/2022/mm12345-1/processing/12345_remap.nxs /dls/i16/data/2022/mm12345-1/12345.nxs
    :param input_files: list of scan file locations
    :param output_directory: str localtion of output files
    :param step: [dh, dk, dl] step size in each direction - size of voxel in reciprocal lattice units
    :return: list of str commands
    """
    input_files = np.asarray(input_files, dtype=str).reshape(-1).tolist()
    return [rsmap_command(file, output_directory, step=step) for file in input_files]


def msmapper_script(input_files, output_file, start=None, shape=None, step=None, bean_file=TEMP_BEAN):
    """
    Create a script that generates a bean file and runs msmapper
     currently only allows a few standard inputs: hkl_start, shape and step values.
    :param input_files: list of scan file locations
    :param output_file: str localtion of output file
    :param start: [h, k, l] start of box
    :param shape: [n, m, o] size of box in voxels
    :param step: [dh, dk, dl] step size in each direction - size of voxel in reciprocal lattice units
    :param bean_file: temporary input file
    :return: str script
    """
    input_files = np.asarray(input_files, dtype=str).reshape(-1).tolist()
    # Remove empty entries
    while '' in input_files:
        input_files.remove('')
    if step is None:
        step = [0.001, 0.001, 0.001]
    else:
        step = np.asarray(step, dtype=float).reshape(-1).tolist()

    if shape is not None:
        shape = np.asarray(shape, dtype=int).reshape(-1).tolist()

    if start is not None:
        start = np.asarray(start, dtype=float).reshape(-1).tolist()

    replace = {
        "inputs": str(list(input_files)),  # Filename of scan file
        "output": "'%s'" % output_file,
        "splitterName": '"gaussian"',
        "splitterParameter": '2.0',
        "scaleFactor": '2.0',
        "step": str(list(step)),
        "start": str(list(start)),
        "shape": str(list(shape)),
        "reduceToNonZero": 'False',
        "bean_file": bean_file,
    }

    with open(TEMPLATE, 'r') as f:
        template = f.read()
    for key, item in replace.items():
        # print("{{%s}}" % key, replace[key], template.count("{{%s}}" % key))
        template = template.replace("{{%s}}" % key, item)
    return template


def plotter_script(output_file):
    """
    Create a script that generates a bean file and runs msmapper
     currently only allows a few standard inputs: hkl_start, shape and step values.
    :param output_file: str localtion of output file
    :return: str script
    """

    replace = {
        "filename": "'%s'" % output_file,
    }

    with open(PLOT_TEMPLATE, 'r') as f:
        template = f.read()
    for key, item in replace.items():
        # print("{{%s}}" % key, replace[key], template.count("{{%s}}" % key))
        template = template.replace("{{%s}}" % key, item)
    return template


def create_bean(input_files, output_file, start=None, shape=None, step=None,
                output_mode=None, to_crystal=None, normalisation=None, polarisation=None,
                detector_region=None, reduce_box=None, third_axis=None,
                azi_plane_normal=None):
    """
    Create a bean file for msmapper in a temporary directory
     currently only allows a few standard inputs: hkl_start, shape and step values.
    :param input_files: list of scan file locations
    :param output_file: str location of output file
    :param start: [h, k, l] start of box (None to omit and calculate autobox)
    :param shape: [n, m, o] size of box in voxels (None to omit and calcualte autobox)
    :param step: [dh, dk, dl] step size in each direction - size of voxel in reciprocal lattice units
    :param output_mode: 'Volume_HKL' or 'Volume_Q' type of calculation
    :param to_crystal: for Volume_Q, use the crystal frame if True, or Lab frame otherwise
    :param normalisation: Monitor value to use for normalisation, e.g. 'rc'
    :param polarisation: Bool apply polarisation correction
    :param detector_region: [sx, ex, sy, ey] region of interest on detector
    :param reduce_box: Bool, reduce box to non-zero elements
    :param third_axis: [h, k, l] direction of Z-axis of voxel grid
    :param azi_plane_normal: [h, k, l] sets X-axis of voxel grid, normal to Z-axis
    :return: str file location of bean file
    """
    input_files = np.asarray(input_files, dtype=str).reshape(-1).tolist()
    # Remove empty entries
    while '' in input_files:
        input_files.remove('')

    version = msmapper_version()

    if step is None:
        step = [0.001, 0.001, 0.001]
    else:
        step = np.asarray(step, dtype=float).reshape(-1).tolist()

    if shape is not None:
        shape = np.asarray(shape, dtype=int).reshape(-1).tolist()

    if start is not None:
        start = np.asarray(start, dtype=float).reshape(-1).tolist()

    if reduce_box is None:
        reduce_box = False

    bean = {
        "inputs": input_files,  # Filename of scan file
        "output": output_file,
        # Output filename - must be in processing directory, or somewhere you can write to
        "splitterName": "gaussian",  # one of the following strings "nearest", "gaussian", "negexp", "inverse"
        "splitterParameter": 2.0,
        # splitter's parameter is distance to half-height of the weight function.
        # If you use None or "" then it is treated as "nearest"
        "scaleFactor": 2.0,
        # the oversampling factor for each image; to ensure that are no gaps in between pixels in mapping
        "step": step,
        # a single value or list if 3 values and determines the lengths of each side of the voxels in the volume
        "start": start,  # location in HKL space of the bottom corner of the array.
        "shape": shape,  # size of the array to create for reciprocal space volume
        "reduceToNonZero": reduce_box  # True/False, if True, attempts to reduce the volume output
    }
    if output_mode:
        bean['outputMode'] = output_mode
    if version > 1.8 and to_crystal is not None:
        bean['toCrystalFrame'] = to_crystal
    if normalisation:
        bean['monitorName'] = normalisation
    if polarisation:
        bean['correctPolarization'] = polarisation
    if detector_region:
        bean['region'] = detector_region
    if version > 1.7 and third_axis:
        bean['thirdAxis'] = np.array(third_axis).tolist()
        bean['aziPlaneNormal'] = np.array(azi_plane_normal).tolist()
    return bean


def create_bean_file(bean_file=None, **kwargs: Unpack[Options]):
    """
    Create a bean file for msmapper in a temporary directory
     currently only allows a few standard inputs: hkl_start, shape and step values.
    :param bean_file: filename of bean file ('bean.json'), or None to create temporary file
    :params: see mapper_runner.create_bean
    :return: str file location of bean file
    """
    bean = create_bean(**kwargs)

    if bean_file is None:
        bean_file = TEMP_BEAN
    json.dump(bean, open(bean_file, 'w'), indent=4)
    print('bean file written to: %s' % bean_file)
    return bean_file


def run_msmapper(**kwargs: Unpack[Options]):
    """
    Create the input file and run msmapper
     currently only allows a few standard inputs: hkl_start, shape and step values.
    :params: see mapper_runner.create_bean
    :return: str file location of bean file
    """
    bean_file = create_bean_file(**kwargs)
    msmapper(bean_file)

