#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Functions for interfacing with the C code.

The C code run_simulation.c is located in the sub-folder c_codes/. 
It performs the core computations for simulating the response of a synapse under a stimulation pattern.
A C code was necessary to speed up the computation time. 

This module uses Python's built-in library ctypes to interface with the C code.
Its uses are the following:
    * importing the C library 
    * getting the objects in the right format for running the C code properly
    * calling the C function and returning the results
"""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"

import os
import ctypes as ct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from subprocess import call
from UBC_synapse import parameters as prm
print("simulate successfully imported")


def compute_coords_measures(S, M_site):
    """Computes the coordinates of the points on the synapse where the glutamate concentrations will be recorded.
        
    Parameters
    ----------
    S : array_like, dtype = int
        Matrix representing the individual synaptic sites areas.
        The synaptic areas are represented by values "1".
        The centers of the synaptic sites remain spotted by values "2".
        The other pixels contain values "0".
    M_site : array-like, dtype = int
        Small patch (matrix) representing a single synaptic site area.
    
    Returns
    ----------
    coords : array-like, dtype = int
        Array of dimension (4, 2), containing the locations of the points where glutamate concentrations are recorded.
        Column 0 : X coordinates.
        Column 1 : Y coordinates.
        Each row corresponds to one point.
        By default, four locations are selected. They are aligned on the middle row of synaptoc sites.
            * 1. point inside a synaptic site, at the center of the synapse
            * 2. point outside a synaptic site, at the center of the synapse
            * 3. point inside a synaptic site, at the periphery of the synapse
            * 4. point outside a synaptic site, at the periphery of the synapse
    """
    # Desired spacing between in and out points:
    spacing = M_site.shape[0]
    # Retrieving the locations of the centers of the synaptic sites:
    sites_coords = np.argwhere(S==2)
    # Retriving the rows of the lines of synaptic centers:
    rows = list(set(site[0] for site in sites_coords))
    rows.sort()
    # Row of the middle line, which crosses the center of the synapse:
    l_center = int(np.median(np.array(rows)))
    # Retrieveing the columns of the centers on this line:
    cols = [site[1] for site in sites_coords if site[0]==l_center]
    # Computing the coordinates:
    c_center_in = int(np.median(np.array(cols)))
    c_center_out = c_center_in + spacing
    c_periph_in = cols[-1]
    c_periph_out = c_periph_in + spacing
    coords = np.concatenate((np.array([l_center]*4).reshape(-1,1), 
                            np.array([c_center_in, c_center_out, c_periph_in, c_periph_out]).reshape(-1,1)),
                            axis=1)
    return coords

def coords_ref(coords) :
    """Definition of reference names for the points where glutamate concentrations are recorded.

    Parameters
    ----------
    coords : array-like, dtype = int
        Array of dimension (4, 2), containing the locations of the points where glutamate concentrations are recorded.        See compute_coords_measures() for details.

    Returns
    ----------
    list of str
        List containing the reference names for the points.
        A reference name is a string that can be used to specify filnames when saving the results.
        Each name is of the shape : 'x-y', with x and y the coordinates of the point.
    """
    return ['{}-{}'.format(coords[i,0], coords[i,1]) for i in range(coords.shape[0])]

def unpack_resglu(resglu, npnts=4):
    """Unpacks the concentrations of glutamate recorded at different points on the synapse.

    Parameters
    ----------
    resglu : array-like, dtype = float
        Array of dimension (npnts*nit,), containing glutamate concentrations. 
        It is part of the raw data returned by running the C code.
        Each bloc of length npnts correspond to one time step, and stores the concentrations at the different locations of recording.
    npnts : int
        Number of points on the synapse where the glutamate concentrations were recorded.

    Returns
    ----------
    glus : list of Series
        List containing the time series of glutamate concentrations at each recording location.
    """
    nit = len(resglu)
    glus = []
    for start in range(npnts):
        glus.append(pd.Series(resglu[np.arange(start,nit, step=npnts)]))
    return glus

# =====================================
# Building an interface with the C code
# =====================================
"""Importation of the C library."""
status = call("gcc -shared -Wl,-soname,run_simulation -o run_simulation.so -fPIC run_simulation.c",cwd="./UBC_synapse/c_code/",shell=True)
run_sim = ct.cdll.LoadLibrary('./UBC_synapse/c_code/run_simulation.so')
run_sim.check_import()
"""Specification of the object types of input and output of the C function run_sim.sim()"""
run_sim.sim.argtypes = [np.ctypeslib.ndpointer(dtype = np.int),
                        np.ctypeslib.ndpointer(dtype = np.int),
                        ct.c_int, 
                        ct.c_int,
                        ct.c_int,
                        np.ctypeslib.ndpointer(dtype = np.int),
                        ct.c_double,
                        ct.c_double,
                        ct.c_double,
                        ct.c_double,
                        np.ctypeslib.ndpointer(dtype = np.int),
                        ct.c_int,
                        ct.c_double,
                        ct.c_double,
                        ct.c_double,
                        ct.c_double,
                        ct.c_double,
                        ct.c_double,
                        np.ctypeslib.ndpointer(dtype = np.double),
                        np.ctypeslib.ndpointer(dtype = np.double),
                        np.ctypeslib.ndpointer(dtype = np.double)]           
run_sim.sim.restype  = ct.c_void_p        
 

def execute_c_code(syn, patt, **kwargs):
    """Calls the C code to run a simulation and prepares the results to be easily used in Python.

    Parameters
    ----------
    syn : Synapse object
        Synapse to be stimulated by the pattern, containing the matrixes among its attributes.
    patt : Pattern object.
        Pattern of stimulation, containing the sequence of neurotransmitter release events among its attributes.
    **kwargs :
        Additional parameters to overwrite the default values.
        coords : array-like, dtype = int
            Array of dimension (4, 2), containing the locations of the points where glutamate concentrations are recorded.
            If not provided, it is automatically computed by compute_coords_measures().
            See compute_coords_measures() for details.
        quant : float
            Quantum of glutamate release (mM) in a spike.
            Default value: 4 mmol (= 4000 umol).
        diff : float
            Diffusion coefficient of glutamate (m2/s).
            Default value: 7.6e-10 m2/s.
        qin : float
            Equivalent of diffusion coefficient adapted in the model's dimensions, in px and tsteps.
            qin applies inside the glomerulus.
            If not provided, qin is computed as diff*tstep/(res**2)
        qout : float
            qout is the same as qin, for outside the synapse.
            If not provided, qout = qin (no diffusion barrier at the limit of glomerulus).

    Returns
    ----------
    glus : list of Series
        List containing the time series of glutamate concentrations at each recording location.
    resAMPAtot_series : Series
        Evolution of the activation state of the AMPA channels (mean activation of the population).
     resV_series : Series
        Evolution of the UBC's membrane potential.
    """
    # Creating arguments from the parameters:
    # I = np.multiply(syn.I, 1).astype(int)
    # S = np.multiply(syn.S, 1).astype(int)
    I = syn.I.copy()
    S = syn.S.copy()
    S[S==2] = 1 # removing the spots of synaptic centers
    nr, nc = I.shape
    spat = patt.spat.copy()
    nit = len(spat)
    tstep = patt.tstep
    coords = kwargs.get('coords', compute_coords_measures(syn.S, syn.M_site))
    npnts = coords.shape[0]
    quant = kwargs.get('quant', 4)
    diff = kwargs.get('diff', 7.6e-10)     
    qin = kwargs.get('qin', diff*patt.tstep/(syn.res**2))
    qout = kwargs.get('qout', qin)
    N_AMPA = syn.N_AMPA
    E_AMPA = syn.E_AMPA
    G_AMPA = syn.G_AMPA
    R = syn.R
    C = syn.C
    V0 = syn.V0
    # Creating variables to store the results generated by the C code:
    resglu = np.zeros(nit*npnts, dtype=np.double)
    resAMPAtot = np.zeros(nit, dtype=np.double)
    resV = np.zeros(nit, dtype=np.double)
    # Running the C code:
    run_sim.sim(I, S, nr, nc,
                nit, spat, tstep, 
                quant, qin, qout, 
                coords, npnts,
                N_AMPA, E_AMPA, G_AMPA, 
                R, C, V0,
                resglu, resAMPAtot, resV)
    # Formatting the data for convenience:
    glus = unpack_resglu(resglu)
    resAMPAtot_series = pd.Series(resAMPAtot)
    resV_series = pd.Series(resV)
    return coords, glus, resAMPAtot_series, resV_series
