# Parametrized Left Ventricle Geometry and Mesh Generator
# Gmsh is used as the backend
#
# Authors
# Bernardo Martins Rocha <bernardomartinsrocha@gmail.com>
# Joventino de Oliveira Campos <joventinoo@gmail.com>
# January, 2019

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def bullseye_plot(ax, data, segBold=None, cmap=None, norm=None,
                  labels=[], labelProps={}):
    """
    Bullseye representation for the left ventricle.

    Parameters
    ----------
    ax : axes
    data : list of int and float
        The intensity values for each of the 17 segments
    segBold: list of int, optional
        A list with the segments to highlight
    cmap : ColorMap or None, optional
        Optional argument to set the desired colormap
    norm : Normalize or None, optional
        Optional argument to normalize data into the [0.0, 1.0] range


    Notes
    -----
    This function create the 17 segment model for the left ventricle according
    to the American Heart Association (AHA) [1]_

    References
    ----------
    .. [1] M. D. Cerqueira, N. J. Weissman, V. Dilsizian, A. K. Jacobs,
        S. Kaul, W. K. Laskey, D. J. Pennell, J. A. Rumberger, T. Ryan,
        and M. S. Verani, "Standardized myocardial segmentation and
        nomenclature for tomographic imaging of the heart",
        Circulation, vol. 105, no. 4, pp. 539-542, 2002.
    """
    if segBold is None:
        segBold = []

    linewidth = 2
    data = np.array(data).ravel()

    if cmap is None:
        cmap = plt.cm.viridis

    if norm is None:
        norm = mpl.colors.Normalize(vmin=data.min(), vmax=data.max())

    theta = np.linspace(0, 2 * np.pi, 768)
    r = np.linspace(0.2, 1, 4)

    # Create the bound for the segment 17
    for i in range(r.shape[0]):
        ax.plot(theta, np.repeat(r[i], theta.shape), '-k', lw=linewidth)

    # Create the bounds for the segments 1-12
    for i in range(6):
        theta_i = np.deg2rad(i * 60)
        ax.plot([theta_i, theta_i], [r[1], 1], '-k', lw=linewidth)

    # Create the bounds for the segments 13-16
    for i in range(4):
        theta_i = np.deg2rad(i * 90 - 45)
        ax.plot([theta_i, theta_i], [r[0], r[1]], '-k', lw=linewidth)

    # Fill the segments 1-6
    r0 = r[2:4]
    r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
    for i in range(6):
        # First segment start at 60 degrees
        theta0 = theta[i * 128:i * 128 + 128] + np.deg2rad(60)
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((128, 2)) * data[i]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if labels:
            ax.annotate(labels[i], xy=(theta0[0,0]+30*np.pi/180,np.mean(r[2:4])),
                        ha='center', va='center', **labelProps)

        if i + 1 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)
            ax.plot(theta0[0], [r[2], r[3]], '-k', lw=linewidth + 1)
            ax.plot(theta0[-1], [r[2], r[3]], '-k', lw=linewidth + 1)

    # Fill the segments 7-12
    r0 = r[1:3]
    r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
    for i in range(6):
        # First segment start at 60 degrees
        theta0 = theta[i * 128:i * 128 + 128] + np.deg2rad(60)
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((128, 2)) * data[i + 6]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if labels:
            ax.annotate(labels[i+6], xy=(theta0[0,0]+30*np.pi/180,np.mean(r[1:3])),
                        ha='center', va='center', **labelProps)

        if i + 7 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)
            ax.plot(theta0[0], [r[1], r[2]], '-k', lw=linewidth + 1)
            ax.plot(theta0[-1], [r[1], r[2]], '-k', lw=linewidth + 1)

    # Fill the segments 13-16
    r0 = r[0:2]
    r0 = np.repeat(r0[:, np.newaxis], 192, axis=1).T
    for i in range(4):
        # First segment start at 45 degrees
        theta0 = theta[i * 192:i * 192 + 192] + np.deg2rad(45)
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((192, 2)) * data[i + 12]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if labels:
            ax.annotate(labels[i+12], xy=(theta0[0,0]+45*np.pi/180,np.mean(r[0:2])),
                        ha='center', va='center', **labelProps)

        if i + 13 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)
            ax.plot(theta0[0], [r[0], r[1]], '-k', lw=linewidth + 1)
            ax.plot(theta0[-1], [r[0], r[1]], '-k', lw=linewidth + 1)

    # Fill the segments 17
    if data.size == 17:
        r0 = np.array([0, r[0]])
        r0 = np.repeat(r0[:, np.newaxis], theta.size, axis=1).T
        theta0 = np.repeat(theta[:, np.newaxis], 2, axis=1)
        z = np.ones((theta.size, 2)) * data[16]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if labels:
            ax.annotate(labels[16], xy=(theta0[0,0]+np.pi/180,0),
                        ha='center', va='center', **labelProps)

        if 17 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)

    ax.set_ylim([0, 1])
    ax.set_yticklabels([])
    ax.set_xticklabels([])

def plot_aha(data):
    
    # Make a figure and axes with dimensions as desired.
    fig, ax = plt.subplots(figsize=(8, 5), nrows=1, ncols=1,
                       subplot_kw=dict(projection='polar'))
    fig.canvas.set_window_title('Left Ventricle WT in AHA Format')

    # Create the axis for the colorbars
    #axl = fig.add_axes([0.14, 0.15, 0.2, 0.05])
    axl = fig.add_axes([0.35, 0.05, 0.35, 0.05])

    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    cmap = mpl.cm.winter
    bounds = np.round(np.linspace(np.min(data), np.max(data),5), decimals=2)
    norm = mpl.colors.Normalize(vmin=np.min(data), vmax=np.max(data), clip=False)
    cmap.set_over(str(np.min(data)))
    cmap.set_under(str(np.max(data)))

    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    cb1 = mpl.colorbar.ColorbarBase(axl, cmap=cmap, norm=norm, orientation='horizontal')
    cb1.set_label("wt")

    # Create the 17 segment model
    lbls =['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
           '11', '12','13', '14', '15', '16', '17']
    bullseye_plot(ax, data, cmap=cmap, norm=norm, labels=lbls)
    ax.set_title('Wall thickness (mm)')

    #fig.savefig('lv_wallthickness.png', inches_bbox='tight')
    plt.show()

# ---------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    
    sample_data = [6.21, 6.38, 6.23, 5.82, 5.38, 6.66,
                   6.40, 7.06, 8.41, 6.64, 5.94, 6.92,
                   5.47, 6.19, 5.40, 5.87,
                   4.37]
    plot_aha(sample_data)



