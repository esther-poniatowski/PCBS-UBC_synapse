#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Functions for generating representations of the synapses' morphology and behaviors."""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"


print("visualize successfully imported")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def print_matrix(S, I, res, print_coords=False, coords=None):
    """Displays the matrixes representing a synapse.
 
    Parameters
    ----------
    S : array_like, dtype = int
        Matrix representing the individual synaptic sites areas.
        The synaptic areas are represented by values "1".
        The centers of the synaptic sites remain spotted by values "2".
        The other pixels contain values "0".
    I : array_like, dtype = int
        Matrix representing the glomerular area (entire synapse).
        It contains values "1" within a disk representing the synapse, and values "0" elsewhere.
    print_coords : bool, optional
        If True, the points where glutamate concentrations were recorded will be displayed.
        Default value: False.
    coords : array-like, dtype = int
        Array of dimension (4, 2), containing the locations of the points where glutamate concentrations are recorded.
        Column 0 : X coordinates.
        Column 1 : Y coordinates.
        Each row corresponds to one point.
    """
    col1 = 'ivory'
    col2 = 'steelblue'
    col3 = 'gold'
    col4 = 'black'
    cmap = [col1, col2, col3]
    patch1 = Patch(color=col2, label='Glomerule')
    patch2 = Patch(color=col3, label='Glutamate release sites')
    patch3 = Patch(color=col1, label='Extraglomerular area')
    patch4 = Patch(color=col4, label='Locations for recording glutamate concentrations')
    handles = [patch1, patch2, patch3]
    extent = [0,S.shape[0]*res*10e3, 0,S.shape[0]*res*10e3]
    S[S==2] = 1
    P = S + I
    if print_coords :
        for i in range(coords.shape[0]):
            P[coords[i,0], coords[i,1]] = 3
        cmap.append(col4)
        handles.append(patch4)
    img = plt.imshow(P,interpolation='nearest', cmap=ListedColormap(cmap), extent=extent)
    plt.xlabel('µm')
    plt.ylabel('µm')
    plt.legend(handles=handles, bbox_to_anchor=(1.04,1), loc="upper left")
    plt.title("Synapse morphology")
    plt.show()

def print_M_site(M_site) :
    """Displays a single synaptic site.

    Parameters
    ----------
    M_site : array-like, dtype = int
        Small patch (matrix) representing a single synaptic site area.
    """
    _, ax = plt.subplots()
    ax.imshow(M_site, extent=(0, M_site.shape[1], M_site.shape[0], 0))
    ax.set_yticks([0,1,2,3,4])
    ax.grid(color='w', linewidth=2)
    ax.set_frame_on(False)
plt.show()

def print_pattern(spat, tstep, nticks=5) :
    """Displays a stimulation pattern.

    Parameters
    ----------
    spat : array_like, dtype = int
        1D boolean array, encoding the release events by values '1'.
    tstep : float
        Temporal resolution used for the pattern instance, in s. Conversion factor, specifying the duration of 1 time step (iteration) in the model.
    nticks : int
        Number of time marks on the x axis.
        Default value: 5.
    """
    fig, ax = plt.subplots()
    ax.plot(spat, label='Spike events')
    xlab = np.round(np.linspace(0, len(spat)*tstep*10e3, nticks), decimals=0)
    xloc = np.linspace(0, len(spat), nticks)
    ax.set_xticks(xloc)
    ax.set_xticklabels(xlab)
    ax.set_xlabel('Time (ms)')
    ax.get_yaxis().set_visible(False)
    ax.legend(loc='upper left', bbox_to_anchor=(0, -0.1))
    plt.show()

