"""
Plot completed msmapper files
"""

import os
import glob
import time
import numpy as np
import matplotlib.pyplot as plt
import hdfmap


DEFAULT_FONT = 'DejaVu Serif'#'Times New Roman'
DEFAULT_FONTSIZE = 14
FIG_SIZE = [12, 8]
FIG_DPI = 80


def set_plot_defaults(rcdefaults=False):
    """
    Set custom matplotlib rcparams, or revert to matplotlib defaults
    These handle the default look of matplotlib plots
    See: https://matplotlib.org/stable/tutorials/introductory/customizing.html#the-default-matplotlibrc-file
    :param rcdefaults: False*/ True, if True, revert to matplotlib defaults
    :return: None
    """
    if rcdefaults:
        print('Return matplotlib rcparams to default settings.')
        plt.rcdefaults()
        return

    plt.rc('figure', figsize=FIG_SIZE, dpi=FIG_DPI, autolayout=False)
    plt.rc('lines', marker='o', color='r', linewidth=2, markersize=6)
    plt.rc('errorbar', capsize=2)
    plt.rc('legend', loc='best', frameon=False, fontsize=DEFAULT_FONTSIZE)
    plt.rc('axes', linewidth=2, titleweight='bold', labelsize='large')
    plt.rc('xtick', labelsize='large')
    plt.rc('ytick', labelsize='large')
    plt.rc('axes.formatter', limits=(-3, 3), offset_threshold=6)
    # default colourmap, see https://matplotlib.org/stable/gallery/color/colormap_reference.html
    plt.rc('image', cmap='viridis')
    # Note font values appear to only be set when plt.show is called
    plt.rc(
        'font',
        family='serif',
        style='normal',
        weight='bold',
        size=DEFAULT_FONTSIZE,
        serif=[DEFAULT_FONT, 'Times', 'DejaVu Serif']
    )
    # plt.rcParams["savefig.directory"] = os.path.dirname(__file__) # Default save directory for figures


def get_remap(nexus_file):
    """
    Returns HdfMap.NexusLoader object
    :param nexus_file:
    :return: scan
    """
    scan = hdfmap.NexusLoader(nexus_file)
    if 'volume' in scan.map:
        print('File contains remapped data')
        return scan
    # Find remapped files
    path, filename = os.path.split(nexus_file)
    name, ext = os.path.splitext(filename)
    print('Finding remapped files in %s' % path)
    for ntries in range(1):
        # remapping may take some time, so keep checking until finished
        remap_files = glob.glob(path + '/processed/*.nxs')
        rmf = [file for file in remap_files if name in file]
        if len(rmf) > 0:
            scan = hdfmap.NexusLoader(rmf[0])
            if 'volume' in scan.map:
                print(f"Remapped file loaded: {scan.filename}")
                return scan
            print(f"Wait for file to finish writing: {scan.filename}")
            time.sleep(10)  # wait for file to finish writing
        else:
            print('Remapped file does not exist yet, try again in 10 s')
            time.sleep(10)
    raise Exception("Scan does not exist or doesn't contain data")


def plot_scan(nexus_file):
    """Default plot scan"""
    scan = hdfmap.NexusLoader(nexus_file)
    xdata, ydata = scan('axes, signal')
    xlabel, ylabel = scan('_axes, _signal')
    plt.figure()
    plt.plot(xdata, ydata, label=nexus_file)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(scan.format('{filename}: {scan_command}'))
    plt.show()


def slider_scan(nexus_file):
    """Plot scan images with slider"""
    nxs_map = hdfmap.create_nexus_map(nexus_file)
    image_path = nxs_map.get_image_path()
    image_len = nxs_map.datasets[image_path].shape[0]
    if not image_path:
        raise Exception('File contains no image data')

    with hdfmap.load_hdf(nexus_file) as hdf:
        image = nxs_map.get_image(hdf, 0)

        fig = plt.figure()
        axes = fig.add_subplot(111)
        pcm = axes.pcolormesh(image)
        axes.invert_yaxis()
        axes.axis('image')
        axes.set_title(nexus_file)

        # colorbar
        fig.colorbar(pcm, ax=axes, label='Pixel intensity', shrink=0.8)

        # Move axes for slider
        bbox = axes.get_position()
        left, bottom, width, height = bbox.x0, bbox.y0, bbox.width, bbox.height
        change_in_height = height * 0.1
        new_position = [left, bottom + 2 * change_in_height, width, height - 2 * change_in_height]
        new_axes_position = [left, bottom, width, change_in_height]

        axes.set_position(new_position, 'original')
        new_axes = axes.figure.add_axes(new_axes_position)

        sldr = plt.Slider(new_axes, label=image_path, valmin=0, valmax=image_len-1,
                          valinit=0, valfmt='%0.0f')

        def update(val):
            """Update function for pilatus image"""
            imgno = int(round(sldr.val))
            im = nxs_map.get_image(hdf, imgno)
            pcm.set_array(im.flatten())
            plt.draw()
            # fig.canvas.draw()

        sldr.on_changed(update)
        plt.show()


def slider_remap(nexus_file):
    """Plot remapped data as slider"""
    nxs_map = hdfmap.create_nexus_map(nexus_file)
    image_path = nxs_map['volume']
    if not image_path:
        raise Exception('File contains no image data')

    with hdfmap.load_hdf(nexus_file) as hdf:
        image = nxs_map.get_data(hdf, image_path, index=(slice(None), slice(None), 0))
        if 'h_axis' in nxs_map:
            x_axis, y_axis, z_axis = nxs_map.eval(hdf, 'h_axis, k_axis, l_axis')
            x_label, y_label, z_label = 'H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)'
        elif 'x_axis' in nxs_map:
            x_axis, y_axis, z_axis = nxs_map.eval(hdf, 'x_axis, y_axis, z_axis')
            x_label, y_label, z_label = u'Qx [\u212B\u207B\u00B9]', u'Qy [\u212B\u207B\u00B9]', u'Qz [\u212B\u207B\u00B9]'
        else:
            raise Exception(f"axes are not recognised in file: {nexus_file}")
        yy, xx = np.meshgrid(y_axis, x_axis)

        fig = plt.figure()
        axes = fig.add_subplot(111)
        pcm = axes.pcolormesh(xx, yy, image)
        axes.set_title(nexus_file)
        axes.axis('image')
        axes.set_xlabel(x_label)
        axes.set_ylabel(y_label)

        # colorbar
        fig.colorbar(pcm, ax=axes, label='Voxel intensity', shrink=0.8)

        # Move axes for slider
        bbox = axes.get_position()
        left, bottom, width, height = bbox.x0, bbox.y0, bbox.width, bbox.height
        change_in_height = height * 0.1
        new_position = [left, bottom + 2 * change_in_height, width, height - 2 * change_in_height]
        new_axes_position = [left, bottom, width, change_in_height]

        axes.set_position(new_position, 'original')
        new_axes = axes.figure.add_axes(new_axes_position)

        sldr = plt.Slider(new_axes, label=z_label, valmin=z_axis[0], valmax=z_axis[-1], valinit=z_axis[0])

        def update(val):
            """Update function for pilatus image"""
            idx = np.argmin(np.abs(sldr.val - z_axis))
            im = nxs_map.get_data(hdf, image_path, index=(slice(None), slice(None), idx))
            pcm.set_array(im.flatten())
            plt.draw()
            # fig.canvas.draw()

        sldr.on_changed(update)
        plt.show()


def plot_remap_hkl(nexus_file):
    """
    Plot Remapped file
    :param nexus_file: str
    :return:
    """

    scan = get_remap(nexus_file)

    # Get reciprocal space data from file
    if 'h_axis' in scan.map:
        x_axis, y_axis, z_axis = scan('h_axis, k_axis, l_axis')
        x_label, y_label, z_label = 'H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)'
    elif 'x_axis' in scan.map:
        x_axis, y_axis, z_axis = scan('x_axis, y_axis, z_axis')
        x_label, y_label, z_label = u'Qx [\u212B\u207B\u00B9]', u'Qy [\u212B\u207B\u00B9]', u'Qz [\u212B\u207B\u00B9]'
    else:
        raise Exception(f"axes are not recognised in file: {nexus_file}")

    vol = scan('volume')

    # print('\n'.join(scan.string(['h-axis', 'k-axis', 'l-axis', 'volume'])))

    # Plot summed cuts
    plt.figure(figsize=(18, 8), dpi=60)
    title = scan.format('{filename}\n{scan_command}')
    plt.suptitle(title, fontsize=18)

    plt.subplot(131)
    plt.plot(x_axis, vol.sum(axis=1).sum(axis=1))
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel('sum axes [1,2]', fontsize=16)

    plt.subplot(132)
    plt.plot(y_axis, vol.sum(axis=0).sum(axis=1))
    plt.xlabel(y_label, fontsize=16)
    plt.ylabel('sum axes [0,2]', fontsize=16)

    plt.subplot(133)
    plt.plot(z_axis, vol.sum(axis=0).sum(axis=0))
    plt.xlabel(z_label, fontsize=16)
    plt.ylabel('sum axes [0,1]', fontsize=16)

    # Plot summed images
    plt.figure(figsize=(18, 8), dpi=60)
    plt.suptitle(title, fontsize=20)
    plt.subplots_adjust(wspace=0.3)

    plt.subplot(131)
    yy, xx = np.meshgrid(y_axis, x_axis)
    plt.pcolormesh(xx, yy, vol.sum(axis=2))
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    # plt.axis('image')
    # plt.colorbar()

    plt.subplot(132)
    zz, xx = np.meshgrid(z_axis, x_axis)
    plt.pcolormesh(xx, zz, vol.sum(axis=1))
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(z_label, fontsize=16)
    # plt.axis('image')
    # plt.colorbar()

    plt.subplot(133)
    zz, yy = np.meshgrid(z_axis, y_axis)
    plt.pcolormesh(yy, zz, vol.sum(axis=0))
    plt.xlabel(y_label, fontsize=16)
    plt.ylabel(z_label, fontsize=16)
    # plt.axis('image')
    plt.colorbar()

    plt.show()


def plot_remap_3dpoints(nexus_file, cmap_name='Greys', cut_ratios=(1e-3, 1e-2, 1e-1)):
    """
    Plot remapped volume in 3D using points
    :param nexus_file: str nexus filename
    :param cmap_name: str matplotlib colormap
    :param cut_ratios: list of cut-ratios, each cut has a different colour and given as ratio of max intensity
    """
    cmap = plt.get_cmap(cmap_name)
    scan = get_remap(nexus_file)
    title = scan.format('{filename}\n{scan_command}')

    # Get reciprocal space data from file
    if 'h_axis' in scan.map:
        x_axis, y_axis, z_axis = scan('h_axis, k_axis, l_axis')
        x_label, y_label, z_label = 'H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)'
    elif 'x_axis' in scan.map:
        x_axis, y_axis, z_axis = scan('x_axis, y_axis, z_axis')
        x_label, y_label, z_label = u'Qx [\u212B\u207B\u00B9]', u'Qy [\u212B\u207B\u00B9]', u'Qz [\u212B\u207B\u00B9]'
    else:
        raise Exception(f"axes are not recognised in file: {nexus_file}")
    yy, xx, zz = np.meshgrid(y_axis, x_axis, z_axis)

    vol = scan('volume')

    alphas = np.linspace(0.1, 1, len(cut_ratios))
    max_volume = vol.max()

    fig = plt.figure(figsize=(12, 10), dpi=100)
    ax = fig.add_subplot(projection='3d')
    for cut, alpha in zip(cut_ratios, alphas):
        filled = vol > (cut * max_volume)
        color = cmap(alpha, alpha=alpha)
        ax.plot(xx[filled].flatten(), yy[filled].flatten(), zz[filled].flatten(), '.', c=color)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_title(title)

    plt.show()


def plot_remap_voxels(nexus_file, cmap_name='Greys', cut_ratios=(1e-3, 1e-2, 1e-1)):
    """
    Plot remapped volume in 3D using points
    :param nexus_file: str nexus filename
    :param cmap_name: str matplotlib colormap
    :param cut_ratios: list of cut-ratios, each cut has a different colour and given as ratio of max intensity
    """
    cmap = plt.get_cmap(cmap_name)
    scan = get_remap(nexus_file)
    title = scan.format('{filename}\n{scan_command}')

    # Get reciprocal space data from file
    if 'h_axis' in scan.map:
        x_axis, y_axis, z_axis = scan('h_axis, k_axis, l_axis')
        x_label, y_label, z_label = 'H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)'
    elif 'x_axis' in scan.map:
        x_axis, y_axis, z_axis = scan('x_axis, y_axis, z_axis')
        x_label, y_label, z_label = u'Qx [\u212B\u207B\u00B9]', u'Qy [\u212B\u207B\u00B9]', u'Qz [\u212B\u207B\u00B9]'
    else:
        raise Exception(f"axes are not recognised in file: {nexus_file}")
    yy, xx, zz = np.meshgrid(y_axis, x_axis, z_axis)

    vol = scan('volume')

    alphas = np.linspace(0.1, 1, len(cut_ratios))
    max_volume = vol.max()

    fig = plt.figure(figsize=(12, 10), dpi=100)
    ax = fig.add_subplot(projection='3d')
    for cut, alpha in zip(cut_ratios, alphas):
        filled = vol > (cut * max_volume)
        nfilled = np.sum(filled)
        if nfilled / vol.size > 0.5:
            print(f"voxel cut={cut}, filled = {nfilled}({nfilled / vol.size:.3%}), skipping...")
            continue
        else:
            print(f"voxel cut={cut}, filled = {nfilled}({nfilled / vol.size:.3%})")
        color = cmap(alpha, alpha=alpha)
        # rebin transparent voxels for speed
        if cut < 0.1:
            ax.voxels(
                xx[::2, ::2, ::2],
                yy[::2, ::2, ::2],
                zz[::2, ::2, ::2],
                filled=filled[::2, ::2, ::2][:-1, :-1, :-1],
                facecolors=color
            )
        else:
            ax.voxels(xx, yy, zz, filled[:-1, :-1, :-1], facecolors=color)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_title(title)

    plt.show()


def plot_histogram(nexus_file, cmap_name='Greys', cut_ratios=(1e-3, 1e-2, 1e-1)):
    """
    Plot remapped volume in 3D using points
    :param nexus_file: str nexus filename
    :param cmap_name: str matplotlib colormap
    :param cut_ratios: list of cut-ratios, each cut has a different colour and given as ratio of max intensity
    """
    cmap = plt.get_cmap(cmap_name)
    scan = get_remap(nexus_file)
    title = scan.format('{filename}\n{scan_command}')
    vol = scan.get_image(index=())
    alphas = np.linspace(0.1, 1, len(cut_ratios))
    max_volume = vol.max()

    ax = plt.figure().add_subplot()
    n, bins, patches = ax.hist(np.log10(vol[vol > 0].flatten()), 100)

    for cut, alpha in zip(cut_ratios, alphas):
        logval = np.log10(cut * max_volume)
        ax.axvline(logval, c='k')
        for patch in patches:
            if patch.xy[0] >= logval:
                patch.set_color(cmap(alpha))

    ax.set_xlabel('Log$_{10}$ Voxel Intensity')
    ax.set_ylabel('N')
    ax.set_title(title)

    plt.show()


def plot_remap_lab(nexus_file, coordinates=None):
    """
    Plot axes of re-mapped nexus file in lab coordinates
    """
    scan = get_remap(nexus_file)
    title = scan.format('{filename}\n{scan_command}')

    wavelength = scan('incident_wavelength') * 10  # incident wavelength in Angstroms
    ubmatrix = scan('ub_matrix') * 2 * np.pi
    corners = -1 * np.array([
        [0, 0, 0], [1, 0, 0], [1, 0, 1], [1, 1, 1],
        [1, 1, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1],
        [1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1],
        [0, 1, 1], [0, 1, 0], [0, 0, 0], [0, 0, 1]
    ])

    # Get reciprocal space data from file
    if 'h_axis' in scan.map:
        h_axis, k_axis, l_axis = scan('h_axis, k_axis, l_axis')
        box = np.array([(h_axis[ii], k_axis[jj], l_axis[kk]) for ii, jj, kk in corners])
        q_space = np.dot(ubmatrix, box.T).T
    elif 'x_axis' in scan.map:
        x_axis, y_axis, z_axis = scan('x_axis, y_axis, z_axis')
        q_space = np.array([(x_axis[ii], y_axis[jj], z_axis[kk]) for ii, jj, kk in corners])
    else:
        raise Exception(f"axes are not recognised in file: {nexus_file}")

    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(projection='3d')

    x_index, y_index, z_index = 0, 2, 1
    labels = [
        'Q || labx',
        'Q || laby',
        'Q || beam'
    ]

    # Beam paths to detector corners
    beam = [0, 0, -2 * np.pi / wavelength]  # in inverse-Angstrom wavevector units
    beam = [0, 0, -1]
    ax.plot([beam[x_index], 0], [beam[y_index], 0], [beam[z_index], 0], 'k-', lw=2)  # incident beam
    for idx in range(4):
        ax.plot([0, q_space[idx, x_index]], [0, q_space[idx, y_index]], [0, q_space[idx, z_index]], 'k-', lw=2)

    # Volume box
    ax.plot(q_space[:, x_index], q_space[:, y_index], q_space[:, z_index], 'k-')

    # crystal axes
    a_axes = np.dot(ubmatrix, [1, 0, 0])
    b_axes = np.dot(ubmatrix, [0, 1, 0])
    c_axes = np.dot(ubmatrix, [0, 0, 1])
    ax.plot([0, a_axes[x_index]], [0, a_axes[y_index]], [0, a_axes[z_index]], 'b-', lw=2)
    ax.plot([0, b_axes[x_index]], [0, b_axes[y_index]], [0, b_axes[z_index]], 'g-', lw=2)
    ax.plot([0, c_axes[x_index]], [0, c_axes[y_index]], [0, c_axes[z_index]], 'y-', lw=2)

    ax.set_xlabel(labels[x_index])
    ax.set_ylabel(labels[y_index])
    ax.set_zlabel(labels[z_index])
    ax.set_title(title)
    ax.set_box_aspect([1, 1, 1])

    if coordinates is not None:
        q_space = np.dot(ubmatrix, coordinates.T).T  # lab frame

        # Volume box
        ax.plot(q_space[:, x_index], q_space[:, y_index], q_space[:, z_index], 'r-')
        ax.plot(q_space[:5, x_index], q_space[:5, y_index], q_space[:5, z_index], 'r-')
        ax.plot(q_space[-5:, x_index], q_space[-5:, y_index], q_space[-5:, z_index], 'r-')

    plt.show()


# Functions from Dans_Diffraction
def bmatrix(a, b=None, c=None, alpha=90., beta=90., gamma=90.):
    """
    Calculate the B matrix as defined in Busing&Levy Acta Cyst. 22, 457 (1967)
    Creates a matrix to transform (hkl) into a cartesian basis:
        (qx,qy,qz)' = B.(h,k,l)'       (where ' indicates a column vector)

    The B matrix is related to the reciprocal basis vectors:
        (astar, bstar, cstar) = 2 * np.pi * B.T
    Where cstar is defined along the z axis

    :param a: lattice parameter a in Anstroms
    :param b: lattice parameter b in Anstroms
    :param c: lattice parameter c in Anstroms
    :param alpha: lattice angle alpha in degrees
    :param beta: lattice angle beta in degrees
    :param gamma: lattice angle gamma in degrees
    :returns: [3x3] array B matrix in inverse-Angstroms (no 2pi)
    """
    if b is None:
        b = a
    if c is None:
        c = a
    alpha1 = np.deg2rad(alpha)
    alpha2 = np.deg2rad(beta)
    alpha3 = np.deg2rad(gamma)

    beta1 = np.arccos((np.cos(alpha2) * np.cos(alpha3) - np.cos(alpha1)) / (np.sin(alpha2) * np.sin(alpha3)))
    beta2 = np.arccos((np.cos(alpha1) * np.cos(alpha3) - np.cos(alpha2)) / (np.sin(alpha1) * np.sin(alpha3)))
    beta3 = np.arccos((np.cos(alpha1) * np.cos(alpha2) - np.cos(alpha3)) / (np.sin(alpha1) * np.sin(alpha2)))

    b1 = 1 / (a * np.sin(alpha2) * np.sin(beta3))
    b2 = 1 / (b * np.sin(alpha3) * np.sin(beta1))
    b3 = 1 / (c * np.sin(alpha1) * np.sin(beta2))

    c1 = b1 * b2 * np.cos(beta3)
    c2 = b1 * b3 * np.cos(beta2)
    c3 = b2 * b3 * np.cos(beta1)

    bm = np.array([
        [b1, b2 * np.cos(beta3), b3 * np.cos(beta2)],
        [0, b2 * np.sin(beta3), -b3 * np.sin(beta2) * np.cos(alpha1)],
        [0, 0, 1 / c]
    ])
    return bm


def latpar2uv(a, b, c, alpha, beta, gamma):
    # From http://pymatgen.org/_modules/pymatgen/core/lattice.html
    alpha_r = np.radians(alpha)
    beta_r = np.radians(gamma)
    gamma_r = np.radians(beta)
    val = (np.cos(alpha_r) * np.cos(beta_r) - np.cos(gamma_r)) \
          / (np.sin(alpha_r) * np.sin(beta_r))
    # Sometimes rounding errors result in values slightly > 1.
    val = abs(val)
    gamma_star = np.arccos(val)
    aa = [a * np.sin(beta_r), a * np.cos(beta_r), 0.0]
    bb = [0.0, b, 0.0]
    cc = [-c * np.sin(alpha_r) * np.cos(gamma_star),
          c * np.cos(alpha_r),
          c * np.sin(alpha_r) * np.sin(gamma_star)]

    return np.round(np.array([aa, bb, cc]), 8)


def RcSp(UV):
    """
    Generate reciprocal cell from real space unit vecors
    Usage:
    UVs = RcSp(UV)
      UV = [[3x3]] matrix of vectors [a,b,c]
    """
    UVs = 2 * np.pi * np.linalg.inv(UV).T
    return UVs


def genq(h, k, l, UVstar):
    avec, bvec, cvec = UVstar
    kk, hh, ll = np.meshgrid(k, h, l)
    shape = np.shape(hh)
    q = avec * hh.reshape(-1, 1) + bvec * kk.reshape(-1, 1) + cvec * ll.reshape(-1, 1)
    qx = q[:, 0].reshape(shape)
    qy = q[:, 1].reshape(shape)
    qz = q[:, 2].reshape(shape)
    return qx, qy, qz


def cal2theta(qmag, energy_kev=17.794):
    """
    Calculate theta at particular energy in keV from |Q|
     twotheta = cal2theta(Qmag,energy_kev=17.794)
    """
    h = 6.62606868E-34  # Js  Plank consant
    c = 299792458  # m/s   Speed of light
    e = 1.6021733E-19  # C  electron charge
    energy = energy_kev * 1000.0  # energy in eV
    # Calculate 2theta angles for x-rays
    twotheta = 2 * np.arcsin(qmag * 1e10 * h * c / (energy * e * 4 * np.pi))
    # return x2T in degrees
    twotheta = twotheta * 180 / np.pi
    return twotheta


def plot_remap_q(nexus_file):
    """
    Plot cuts in Q
    :param nexus_file:
    :return:
    """

    scan = get_remap(nexus_file)

    # Get reciprocal space data from file
    vol = scan('volume')
    if 'h_axis' in scan.map:
        h_axis, k_axis, l_axis = scan('h_axis, k_axis, l_axis')
        # Convert to Q
        a, b, c, alpha, beta, gamma = scan('unit_cell')
        astar, bstar, cstar = RcSp(latpar2uv(a, b, c, alpha, beta, gamma))
        qx, qy, qz = genq(h_axis, k_axis, l_axis, [astar, bstar, cstar])
    elif 'x_axis' in scan.map:
        qx, qy, qz = scan('x_axis, y_axis, z_axis')
    else:
        raise Exception(f"axes are not recognised in file: {nexus_file}")

    # Plot summed images
    plt.figure(figsize=[18, 8], dpi=60)
    title = scan.format('{filename}\n{scan_command}')
    plt.suptitle(title, fontsize=20)
    plt.subplots_adjust(wspace=0.3)

    plt.subplot(131)
    plt.pcolormesh(qx.mean(axis=2), qy.mean(axis=2), vol.sum(axis=2))
    plt.xlabel('Qx [A$^{-1}$]', fontsize=16)
    plt.ylabel('Qy [A$^{-1}$]', fontsize=16)
    # plt.axis('image')
    # plt.colorbar()

    plt.subplot(132)
    plt.pcolormesh(qx.mean(axis=1), qz.mean(axis=1), vol.sum(axis=1))
    plt.xlabel('Qx [A$^{-1}$]', fontsize=16)
    plt.ylabel('Qz [A$^{-1}$]', fontsize=16)
    # plt.axis('image')
    # plt.colorbar()

    plt.subplot(133)
    plt.pcolormesh(qy.mean(axis=0), qz.mean(axis=0), vol.sum(axis=0))
    plt.xlabel('Qy [A$^{-1}$]', fontsize=16)
    plt.ylabel('Qz [A$^{-1}$]', fontsize=16)
    # plt.axis('image')
    plt.colorbar()
    plt.show()


def generate_twotheta_profile(nexus_file):
    """
    Generate a 1D profile of the reciprocal space map in units of wavevector-transfer (|Q|) and angle (two-theta)
    :param nexus_file: filename of the re-mapped nexus file
    :returns: array(two-theta), array(Q), array(intensity)
    """

    scan = get_remap(nexus_file)

    # Get reciprocal space data from file
    energy = scan('incident_energy')
    vol = scan('volume')
    if 'h_axis' in scan.map:
        h_axis, k_axis, l_axis = scan('h_axis, k_axis, l_axis')
        # Convert to Q
        a, b, c, alpha, beta, gamma = scan('unit_cell')
        astar, bstar, cstar = RcSp(latpar2uv(a, b, c, alpha, beta, gamma))
        qx, qy, qz = genq(h_axis, k_axis, l_axis, [astar, bstar, cstar])
    elif 'x_axis' in scan.map:
        qx, qy, qz = scan('x_axis, y_axis, z_axis')
    else:
        raise Exception(f"axes are not recognised in file: {nexus_file}")

    qmag = np.sqrt(qx ** 2 + qy ** 2 + qz ** 2)
    qmag = qmag.reshape(-1)
    qvol = vol.reshape(-1)
    bin_cen = np.arange(qmag.min(), qmag.max(), 0.001)
    bin_edge = bin_cen + 0.005
    bin_pos = np.digitize(qmag, bin_edge) - 1
    bin_sum = [np.mean(qvol[bin_pos == n]) for n in range(len(bin_cen))]
    tth = cal2theta(bin_cen, energy)
    return tth, bin_cen, bin_sum


def plot_qmag(nexus_file):
    """
    Plot two0-theta intensity and intensity vs |Q|
    :param nexus_file:
    :return:
    """

    tth, bin_cen, bin_sum = generate_twotheta_profile(nexus_file)
    scan = get_remap(nexus_file)
    title2 = scan.format('{filename}\n{scan_command}')

    plt.figure(dpi=100)
    plt.plot(bin_cen, bin_sum)
    plt.title(title2, fontsize=20)
    plt.xlabel('|Q| A$^{-1}$', fontsize=20)
    plt.ylabel('Intensity', fontsize=20)

    plt.figure(dpi=100)
    plt.plot(tth, bin_sum)
    plt.title(title2, fontsize=20)
    plt.xlabel('Two-Theta [Deg]', fontsize=20)
    plt.ylabel('Intensity', fontsize=20)
    plt.show()

