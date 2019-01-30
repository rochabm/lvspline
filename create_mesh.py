# Parametrized Left Ventricle Geometry and Mesh Generator
# Gmsh is used as the backend
#
# Authors
# Bernardo Martins Rocha <bernardomartinsrocha@gmail.com>
# Joventino de Oliveira Campos <joventinoo@gmail.com>
# January, 2019

import numpy as np
import lvspline as lvs
from plotwt import plot_aha

if __name__ == "__main__":

    # Example 1 - wall thickness (in meters) from Bai et al. ATLAS paper
    lv_ex1 = [0.00621, 0.00638, 0.00623, 0.00582, 0.00538, 0.00666,
              0.00640, 0.00706, 0.00841, 0.00664, 0.00594, 0.00692,
              0.00547, 0.00619, 0.00540, 0.00587,
              0.00437]

    # Example 2 - wall thickness (in meters)
    lv_ex2 = [0.0058, 0.0064, 0.0102, 0.0108, 0.0084, 0.0061,
              0.0092, 0.0081, 0.0102, 0.0108, 0.0082, 0.0081,
              0.0157, 0.0157, 0.0102, 0.0127,
              0.0157]

    # Short and long axes
    short_axis = 0.04
    long_axis = 0.06

    # Changes the WT of AHA-segments 3 and 9 for testing
    #lv_ex1[8] *= 1.5
    #lv_ex1[2] = lv_ex1[8]

    # Plot diagram with wall thickness
    bp = np.array(lv_ex1)
    bp = 1000 * bp[:]
    plot_aha(bp)

    # Creates the geometry and mesh
    lvs.generateLVMesh(short_axis, long_axis, lv_ex1)

