#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Functions for building the morphology of a synapse.

This modules meets two goals:
    1. Determining relevant ranges of parameters to model synapses. They should be close to those observed in vivo.
    2. Building matrixes to represent the morphology of the synapses.

The model of a synapse consists in three matrixes:
    * the synapse area (glomerulus) is a disk,
    * the glutamate release sites are homogenousely distributed on ths area, according to a hexagonal grid,
    * the extra glomerular area, between the borders of the synapse and the limits of the working area.
These matrixes are boolean matrixes, containing 1s in the pixels representing the corresponding areas, 0s elsewhere.
"""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"

import numpy as np
import matplotlib.pyplot as plt

def convert_dist_to_px(dist, res=0.2e-6):
    """Conversion of a distance in meters to a number of pixels in the model.

    Parameters
    ----------
    dist : float
        Distance to be converted, in m.
    res : float, optional
        Sptatial resolution of the model in m. Conversion factor. 
        Default: 0.2e-6 m, corresponding to 1 px <=> 0.2 um.

    Returns
    -------
    int
        Number of pixels corresponding to the given distance in the model (inferior integrer part).
    """
    return int(dist/res)

def convert_area_to_px(area, res=0.2e-6):
    """Conversion of an area in meters to a number of pixels in the model.

    Parameters
    ----------
    area : float
        Area to be converted, in m².
    res : float, optional
        Sptatial resolution of the model, in m. The conversion factor for areas is thus res**2. 
        Default: 0.2e-6 m, corresponding to 1 px <=> 0.2 um.

    Returns
    -------
    int
        Number of pixels corresponding to the given area in the model (inferior integrer part).
    """
    return int(area/res**2)

def compute_n_px_site(n_sites, area_tot_sites, res):
    """Proposes model parameters to represent the synapse, matching as much as possible the constraints provided by the user.
    
    Parameters
    ----------
    n_sites : int
        Desired number of distinct presynaptic glutamate release sites.
    area_tot_sites : float
        Desired cumulative area of glutamate release.
    res : float
        Sptatial resolution of the model in m. Conversion factor. 

    Returns
    -------
    n_px_site : int
        Number of pixels computed for an individual glutamate release site.
    area_updated : float
        Cumulative area of glutamate release in meters, in a model featuring the proposed parameters.
        It may differ from the target area provided by the user, because of the constraint of working with pixels.
    
    Raises
    ------
    TypeError
        If the number of sites provided is not an integrer.
    """
    assert isinstance(n_sites, int), "n_sites should be an integrer."
    area_px = convert_area_to_px(area_tot_sites, res)
    n_px_site = int(np.ceil(area_px / n_sites))
    area_updated = (n_px_site*n_sites)*(res**2)
    return n_px_site, area_updated

def compute_patch_side(n_px_site):
    """Computes the dimension in pixels for a small patch (matrix) representing a single glutamate release site.

    For convenience, the patches matrixes are squared matrixes.
    When the number of pixels provided is a square number, the computation of the dimension is exact. The future matrix will contain only 1s.
    When the number of pixels provided is not a square number, the function takes the superior integrer part. The future matrix will contain 1s and 0s (in excess pixels).

    Parameters
    ----------
    n_px_site : int
        Number of pixels in a single glutamate release site.

    Returns
    -------
    int
        Dimension in pixels for the matrix representing a glutamate release site (number of pixels in a side).

    See also
    --------
    create_site - Function buidling the corresponding small patch.
    """
    return int(np.ceil(n_px_site**0.5))

def create_site(n_px_site):
    """Creates the small patch (matrix) representing a single synaptic site area.

    To ensure homogeneity between models featuring different individual synaptic areas, a synaptic site matrix is built iteratively.
    The matrix is initialized with a value "1" at the center. At each iteration, a value "1" is added in a new pixel.
    The building method consists in a spiral growth in the trigonometric direction. The order of the moves is Noth, West, South, East.
    At the end of the process, the center remains spotted by a value "2".

    Parameters
    ----------
    n_px_site : int
        Number of pixels in a single glutamate release site.
    
    Returns
    -------
    M_site : array-like, dtype = int
        Small patch (matrix) representing a single synaptic site area.
    """
    c = compute_patch_side(n_px_site)
    M_site = np.zeros((c, c), dtype=int)
    center = (int(c/2), int(c/2))
    M_site[center] = 1
    loc = center
    # Sequence of coordinates moves to make a spiral growth (N:(1,0), W:(0,-1), S:(-1,0), E:(0,1)):
    moves = [([(1,0)]*i + [(0,1)]*i)*(i%2==0) + ([(-1,0)]*i + [(0,-1)]*i)*(i%2==1) for i in range(1, n_px_site+1)]
    # Unpacking to get tuples:
    moves = [x for elem in moves for x in elem]
    for i in range(n_px_site):
        M_site[loc] += 1
        loc = (loc[0]+moves[i][0], loc[1]+moves[i][1])
    return M_site

def compute_dist_btw_sites(d):
    """Computes a reference distance in pixels between lines of synaptic sites centers on the future hexagonal grid.

    For convenience in working on matrixes, the distribution of the individual synaptic sites on the entire synapse is described by two parameters: 
        * the number of columns between two centers along a line,
        * the number of rows separating two *lines* of staggered centers.
    Two constraints bear upon the choice of these parameters:
        * in a hexagonal grid, all synaptic centers must be at equal distance from eachother,
        * the synaptic sites areas should not overlap.
    Pythagore' theorem shows that these constraints impose the number of columns. It should be at least equal to the side of the small patch (matrix) representing a single synaptic site.
    This function computes the adapted number of rows corresponding to the number of columns.
    
    Parameters
    ----------
    d : int
        Number of columns separating two synaptic sites centers along a line.
    
    Returns
    -------
    h : int
        Number of rows separating two lines of synaptic centers. 
    """
    h = int(np.floor((d**2 + (d/2)**2)**0.5))
    return h

def compute_hexgrid_params(n_sites):
    """Computes the parameters of the heganonal grid, to be closer as possible to the desired number of individual synaptic sites.

    A hexagonal grid is a symetric network of points, which are all at equal distance from eachother.
    To gain flexibility in the different number of points that can be obtained, several geometries are consideredfor the grid:
        * "standard hexagonal grid": grid featuring a hexagonal external shape, i.e. whose borders form a perfect hexagon.
        * "derivatives" are obtained by adding points around each side of the hexagon, ensuring symetry:
            * "star shaped grid": obtained from a "standard grid" whose number of points on a side is even, by adding one single point at a distance from the middle of each side ; 
            * "flower shaped grid":
                * either obtained from a "standard grid" whose number of points on a side is odd, 
                * or obtained from a "star shaped grid", by adding two points parallelly to the hexagon's side ;
                * or obtained from another "flower shaped grid", by expanding the "petals" width towards an hexagonal shape.
    In this module, a hexagonal grid is described by two parameters:
        * order: maximal number of hexagons (centered on the grid center) contained in the grid.
            It is also equal to the number of points forming a side of a "standard hexagonal gird".
            Examples : 
            The order 1 encompasses the "standard hexagonal grid" with 7 points, and the "star shaped" grid with 13 points.
            The order 2 encompasses the "standard hexagonal grid" with 19 points, and the "flower shaped" grid with 31 points.
        * p : number of points added to the "standard grid" of the current order, to obtain a particular derivative grid.
            All standard hexagonal grids are associated with p = 0.
            All "star shaped" grids are associated with p = 1.
            The "flower shaped" grids are associated with p >= 1 when their order is odd, p >=2 when their order is even.

    Parameters
    ----------
    n_sites : int
        Number of individual glutamate release sites.
    
    Returns
    -------
    n_points : int
        Proposed number of points, that best fit the desired number of points, while allowing to build a symetric grid.
    order : int
        Oder of the target grid.
        Maximal number of hexagons (centered on the grid center) contained in the grid.
        See above for examples.
    p_max : int
        Reference for the shape of the target grid. 
        Number of points added to the "standard hexagonal grid" of the same order, to obtain the target grid.
        See above for more precise definitions of the geometries.
    """
    n_points = 1
    order = 1
    p_max = 0
    while n_points < n_sites:
        order += 1
        n_points += 6
        p_sup = int(np.floor(order/2))
        p = 0
        if (n_points < n_sites and order % 2 == 0):
            p = 1
            n_points += 6
        while (n_points < n_sites and p < p_sup):
            p += 1
            n_points += 2*6
        p_max = p
    return n_points, order, p_max

def compute_hexgrid_coords(order, p, d):
    """Computes the coordinates of the individual synaptic sites on the hexagonal grid.

    The method for building a grid procedes in several steps:
        1. Alternating two sequences of x coordinates along the y axis.
            The two x sequences differ by an offset which corresponds to half the distance between sites' centers along a line.
            They are separated in the y direction by a margin which corresponds to the distance between *lines* of staggered centers.
        2. Removing external points on each line to match the desired geometry.

    Parameters
    ----------
    order : int
        Oder of the target grid.
        Maximal number of hexagons (centered on the grid center) contained in the grid.
    p : int
        Reference for the shape of the target grid. 
        Number of points added to the standard hexagonal grid of the same order, to obtain the target grid.
    d : int
        Number of columns separating two synaptic sites centers along a line.

    Returns
    -------
    coords_list : list of tuples of int
        Coordinates of the synaptic centers on the hexagonal grid. 
        The center of the grid is assumed to be (0, 0), thus the coordinates can be positive or negative.

    See also
    --------
    compute_hexgrid_params - For a more precise definition of the arguments "order" and "p".
    compute_dist_btw_sites - For a more precise description of the arguments "d".
    """
    h = compute_dist_btw_sites(d)
    # Sequences of x coordinates, deiffierent by an offset of d/2:
    X1 = np.arange(0, order*d, step=d)
    X2 = X1[:len(X1)] + int(np.ceil(d/2))
    # Alernation of the two previous sequences:
    X_list = [X1 if l%2==0 else X2 for l in range(order+1)]
    # Number of derivatives grids which can be obtained from the standard hexagonal grid of the current order:
    p_sup = int(np.floor(order/2))
    # Computation of the threshold, for each line, for which external points have to be removed from the x coordinate sequence: 
    if order%2 == 0:
        p_thres = np.array([i for i in reversed(range(1,p_sup+2))] + [i for i in range(2,p_sup+2)])
    else :
        p_thres = np.array([i for i in reversed(range(1,p_sup+2))] + [i for i in range(1,p_sup+2)])
    # Number of points to be removed on each line: 
    rem_number = np.array([0] + [l-(l%2==0) for l in range(1,order)] + [order-p]) 
    rem_number = rem_number - (p_thres<=p)
    # Building the coordinates line by line:
    coords_list = []
    for l in range(order+1):
        y = l*h
        X = list(X_list[l])
        if rem_number[l] != 0:
            for i in range(rem_number[l]) :
                X.pop()
        for x in X :
            # Generation of four points, to ensure radial symetry around the grid center:
            coords_list.extend([(-x, -y), (-x, y), (x, -y), (x, y)])
    # Removal of double coordinates:
    coords_list = list(dict.fromkeys(coords_list))
    return coords_list

def compute_dimensions(coords_list, c, border_intra, border_extra):
    """Computes the dimension of the final working matrix and the radius of the glomerulus (entire synapse).
    
    The constraints are the following:
        * The glomerulus should contain all the synaptic sites.
        * The glomerulus should include a margin between the most remote synaptic site and the border of the glomerulus.
        * The final working matrix should include a margin of extra-glomerular area, between the border of the glomerulus and the limit of the matrix.
    
    Parameters
    ----------
    coords_list : list of tuples of int
        Coordinates of the synaptic centers on an hexagonal grid. 
        The center of the grid is assumed to be (0, 0), thus the coordinates can be positive or negative.
    c : int
        Dimension of the small patch (matrix) representing a single synaptic site area.
    border_intra : int
        Margin in px between the most external synaptic area and the synapse border.
    border_extra : int
        Margin in px between the synapse border and the limit of the working matrix.

    Returns
    -------
    half_dim : int
        Half dimension of the final matrix.
        More precisely, the number of pixels between the center and the border of the matrix.
        The dimension of the matrix will be 2*half_dim + 1.
    r : int 
        Radius of the glomerulus (entire synapse) in px.

    """
    # Remoteness of the most external point:
    dist_max = max([int(np.ceil((coord[0]**2 + coord[1]**2)**0.5)) for coord in coords_list])
    # Adding the half diagonal of a synaptic site matrix, and the internal margin of the glomerulus:
    r = dist_max + int((c/2**0.5)) + border_intra
    half_dim = r + border_extra
    return half_dim, r

def create_S(coords_list, M_site, border_intra, border_extra):
    """Creates the final matrix representing the individual synaptic sites areas.

    The process is achieved in several steps :
        1. Generating a null matrix of the appropriate dimension.
        2. Translating the coordinated of the synaptic centers, so that the point (0, 0) coincides with the center of the working matrix.
        3. On each location of a synaptic center, "pasting" the small patch (matrix) representing single synaptic site with values "1".

    Parameters
    ----------
    coords_list : list of tuples of int
        Coordinates of the synaptic centers on an hexagonal grid. 
        The center of the grid is assumed to be (0, 0), thus the coordinates can be positive or negative.
    M_site : array-like, dtype = int
        Small patch (matrix) representing a single synaptic site area.
    border_intra : int
        Margin in px between the most external synaptic area and the synapse border.
    border_extra : int
        Margin in px between the synapse border and the limit of the working matrix.

    Returns
    -------
    S : array_like, dtype = int
        Final matrix representing the individual synaptic sites areas.
        The synaptic areas are represented by values "1".
        The centers of the synaptic sites remain spotted by values "2".
        The other pixels contain values "0".
    """
    half_dim = compute_dimensions(coords_list, M_site.shape[0], border_intra, border_extra)[0]
    S = np.zeros((2*half_dim + 1, 2*half_dim + 1), dtype=int)
    # Translation of the coordinates in the positive range:
    for coord in coords_list :
        S[coord[0]+half_dim, coord[1]+half_dim] += 1
    # Distances around each synaptic center, necessary to contain the patch matrix:
    center = np.argwhere(M_site==2)[0]
    xleft = center[0]
    xright = M_site.shape[0] - center[0]
    yleft = center[1]
    yright = M_site.shape[1] - center[1]
    # Retriving the synaptic centers and pastiing the patches:
    anchors = np.argwhere(S == 1)
    for a in anchors :
        S[a[0]-xleft:a[0]+xright, a[1]-yleft:a[1]+yright] = M_site
    return S

def create_I_O(half_dim, r):
    """Creates the final matrixes representing the glomerulus area (entire synapse) and the extraglomerulat working area.

    Parameters
    ----------
    half_dim : int
        Half dimension of the final matrix.
        More precisely, the number of pixels between the center and the border of the matrix.
        The dimension of the matrix will be 2*half_dim + 1.
    r : int 
        Radius of the glomerulus (entire synapse) in px.

    Returns
    -------
    I : array_like, dtype = int
        Matrix representing the glomerular area (entire synapse).
        It contains values "1" within a disk representing the synapse, and values "0" elsewhere.
    O : array_like, dtype = int
        Matrix representing the etraglomerular area. O is the mirror of I, obtained by 1 - I.
    """
    x = np.arange(-half_dim, half_dim+1)
    xv, yv = np.meshgrid(x, x)
    D = (xv**2 + yv**2)**0.5
    I = (D<=r).astype(int)
    O = np.ones((2*half_dim + 1, 2*half_dim + 1), dtype=int) - I
    return I, O

def compute_scaled_d(c, dim, dim_ref):
    """Adapts the distance between synaptic centers alon g aline, in order to get a matrix of a desired dimension.

    This function implements a scaling with a cross products.
    
    Parameters
    ----------
    c : int
        Dimension of a small patch (matrix) representing an individual synaptic site (number of pixels in a side).
    dim : int
        Target dimension of the desired matrix.
    dim_ref : int
        Dimension of a reference matrix, previousely computed to satisfy the geometric contraints imposed by the dimension of the patch matrix.

    Returns
    -------
    int
        Number of columns between synaptic centers along a line in the target matrix.
    """
    return int(np.ceil(c*dim/dim_ref))


def create_all_matrixes(**kwargs):
    """Creates the three final matrixes representing a synapse.

    Parameters
    ----------
    **kwargs : 
        Any parameter to constrain the model.
        dim : int
            Desired diameter for the target glomerulus, in px.
            Default value: None.
            In this case, a reference dimension is computed.
         n_sites : int
            Number of individual glutamate release sites.
            Default value: 37.
        area_tot_sites : float
            Desired cumulative area of glutamate release, in m².
            Default value: 25e-12 m² (between 12 and 40 um², in vivo references given by Mungi).
        res : float
            Sptatial resolution of the model in m. Conversion factor. 
            Default value: 0.2e-6 m, corresponding to 1 px <=> 0.2 um.
        border_intra : int
            Margin in px between the most external synaptic area and the synapse border.
            Default value: 2 px.
        border_extra : int
            Margin in px between the synapse border and the limit of the working matrix.
            Default value: 2 px.
        verbose : bool, optional
            If True, the updated/suggested parameters and information will be printed.
            Default value: False.

    Returns
    -------
    S : array_like, dtype = int
        Matrix representing the individual synaptic sites areas.
        The synaptic areas are represented by values "1".
        The centers of the synaptic sites remain spotted by values "2".
        The other pixels contain values "0".
    I : array_like, dtype = int
        Matrix representing the glomerular area (entire synapse).
        It contains values "1" within a disk representing the synapse, and values "0" elsewhere.
    O : array_like, dtype = int
        Matrix representing the extraglomerulat area. O is the mirror of I, obtained by 1 - I.
    M_site : array-like, dtype = int
        Small patch (matrix) representing a single synaptic site area.
    n_sites_updated : int
        Retained number of pixels computed for an individual glutamate release site.
        Note: It may differ from the target number provided by the user.
    area_updated : float
        Retained cumulative area of glutamate release, in m, with the proposed parameters computed to fit several constraints.
        Note: It may differ from the target area provided by the user.
    dim_updated : int
        Retained dimension of the glumerile, in pixels.
        Note: It may differ from the target area provided by the user.
    """
    # Retriving the specified parameters, setting the default values otherwise:
    dim = kwargs.get("dim", None)
    n_sites = kwargs.get("n_sites", 37)
    area_tot_sites = kwargs.get('area_tot_sites', 25e-12)
    res = kwargs.get('res', 0.2e-6)
    border_intra = kwargs.get('border_intra', 2)
    border_extra = kwargs.get('border_extra', 2)
    # Computation of an appropriate number of sites:
    n_sites_updated, order, p = compute_hexgrid_params(n_sites)
    n_px_site, area_updated = compute_n_px_site(n_sites_updated, area_tot_sites, res)
    c = compute_patch_side(n_px_site)
    # Computation of a reference dimension for the matrix, to contain all the sites:
    # contraint : d = c
    coords_list = compute_hexgrid_coords(order, p, c)
    half_dim, r = compute_dimensions(coords_list, c, border_intra, border_extra)
    dim_updated = 2*r + 1
    # Computation of the final parameters to match the desired dimension as much as possible (if provided):
    if dim != None and dim != dim_updated:
        d = compute_scaled_d(c, dim, dim_updated)
        coords_list = compute_hexgrid_coords(order, p, d)
        half_dim, r = compute_dimensions(coords_list, c, border_intra, border_extra)
        dim_updated = 2*r + 1
    # Creation of the matrixes:
    M_site = create_site(n_px_site)
    S = create_S(coords_list, M_site, border_intra, border_extra)
    I, O = create_I_O(half_dim, r)
    # Printing the results of the computations:
    if kwargs.get("verbose", False) :
        if n_sites_updated != n_sites :
            print("For bluiding a hexgrid, the closest possible value for n_sites is {}.".format(n_sites_updated))
        else :
            print("Number of release sites: {}.".format(n_sites_updated))
        if area_updated != area_tot_sites :
            print("Release area updated with this parameter: {} m^2.".format(area_updated))
        else :
            print("Release area: {} m^2.".format(area_updated))
        print("Area of a single site: {} m^2 ({} px).".format(n_px_site*res, n_px_site))
        print("Diameter obtained for the glomerulus: {} m = {} um ({} px).".format(dim_updated*res, (dim_updated*res)*10e6, dim_updated))
    return S, I, O, M_site, n_sites_updated, area_updated, dim_updated