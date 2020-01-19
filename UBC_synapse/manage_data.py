#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Functions saving and retriving the results."""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"

print("manage_data successfully imported")

import numpy as np
import pandas as pd
import os
import sys
from types import FunctionType
from inspect import getmembers


def path_directory(option, syn=None):
    """Defines a path to a data directory.

    Parameters
    ----------
    option : {'all_syn', 'syn', 'resp'}
        * 'all_syn': to get the path to the top directory 'UBC_data/',
        * 'syn': to get the path to the directory of a synapse 'UBC_data/syn{index}',
        * 'resp': to get the path to the directory of a response of the synapse, stored in the attribute "resp".
    syn : Synapse instance, optional
        Object of the class Synapse, containing the required indexes in its attribute indexes['syn'] and indexes['resp'].
        Necessary if the options 'syn' or 'resp' are used.
        Default value: None.

    Returns
    -------
    path_dir : str
        Absolute path to the desired directory.
    """
    # Initialize indexes if no instance has been saved before:
    if option == 'all_syn':
        extension = ''
    else:
        if syn.indexes['syn'] == None:
            s_ind = 0
        else:
            s_ind = syn.indexes['syn']
        if syn.indexes['resp'] == None:
            r_ind = 0
        else:
            r_ind = syn.indexes['resp']
        extension = {'syn': 'syn{}'.format(s_ind),
                    'resp': os.path.join('syn{}'.format(s_ind), 'resp{}_syn{}'.format(r_ind, s_ind))}[option]
    path_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'UBC_data', extension)
    return path_dir

def path_register(option, syn=None):
    """Defines a path to a register of saves.

    Parameters
    ----------
    option : {'all_syn', 'syn'}
        * 'all_syn': to get the path to the register of all synapses, 'synapses_register.csv',
        * 'syn': to get the path to the register of a synapse 'responses_syn{index}_register.csv'.
    syn : Synapse instance, optional
        Object of the class Synapse, containing the required index in its attribute indexes['syn'].
        Necessary if the options 'syn' or 'resp' are used.
        Default value: None.

    Returns
    -------
    path_dir : str
        Absolute path to the desired register.
    """
    if option == 'all_syn':
        extension = 'synapses_register.csv'
    else:
        extension = 'responses_syn{}_register.csv'.format(syn.indexes['syn'])
    path_reg = os.path.join(path_directory(option, syn), extension)
    return path_reg

def create_directory(option, syn=None):
    """Creates a directory for saving the results, if it does not exist.

    Parameters
    ----------
    option : {'all_syn', 'syn', 'resp'}
        * 'all_syn': to create the top directory 'UBC_data/',
        * 'syn': to create the directory of a synapse 'UBC_data/syn{index}',
        * 'resp': to create the directory of a response of the synapse, stored in the attribute "resp".
    syn : Synapse instance, optional
        Object of the class Synapse, containing the required indexes in its attribute indexes['syn'] and indexes['resp'].
        Necessary if the options 'syn' or 'resp' are used.
        Default value: None.
    """
    try :
        os.makedirs(path_directory(option, syn))
        print("Creating a directory at {}".format(path_directory(option, syn)))
    except FileExistsError:
        print("A directory already exists at {}".format(path_directory(option, syn)))

def create_register(option, syn=None):
    """Creates a register to classify the future saves.

    Parameters
    ----------
    option : {'all_syn', 'syn'}
        * 'all_syn': to create the register of all synapses, 'synapses_register.csv',
        * 'syn': to create the register of the synapse, 'responses_syn{index}_register.csv'.
    syn : Synapse instance, optional
        Object of the class Synapse, containing the required index in its attribute indexes['syn'].
        Necessary if the options 'syn' or 'resp' are used.
        Default value: None.
    """
    columns_opt = {'all_syn': ['dim', 'n_sites', 'area_tot_sites', 
                            'res', 'border_intra', 'border_extra',
                            'N_AMPA', 'G_AMPA', 'E_AMPA', 'V0', 'R', 'C'],
                    'syn': ['mode','tstep', 'dur', 'nit', 
                            'start', 'n_stim', 'f_stim', 
                            'fenv', 'fmin', 'fmax']}
    if not os.path.exists(path_directory(option, syn)): # if directory never created
        create_directory(option, syn)
    if not os.path.exists(path_register(option, syn)): # to avoid overwriting an existing register
        register = pd.DataFrame(columns=columns_opt[option])
        register.index.name = 'Index'
        register.to_csv(path_register(option, syn), float_format='%.16g')
        print("Creating a register at {}.".format(path_register(option, syn)))

def load_register(option, syn=None):
    """Loads a register classifying the previous saves.

    Parameters
    ----------
    option : {'all_syn', 'syn'}
        * 'all_syn': to load the register of all synapses, 'synapses_register.csv',
        * 'syn': to load the register of the synapse, 'responses_syn{index}_register.csv'.
    syn : Synapse instance, optional
        Object of the class Synapse, containing the required index in its attribute indexes['syn'].
        Necessary if the options 'syn' or 'resp' are used.
        Default value: None.

    Returns
    -------
    DataFrame
        If a register exists.
    None
        If no register has been created yet.
    """
    if os.path.exists(path_register(option, syn)) :
        return pd.read_csv(path_register(option, syn), float_precision='high')
    else:
        return None

def filter_register(option, syn=None, **kwargs):
    """Scans a register to identify saves which feature specific properties.

    Parameters
    ----------
    option : {'all_syn', 'syn'}
        * 'all_syn': to scan the register of all synapses, 'synapses_register.csv',
        * 'syn': to scan the register of the synapse, 'responses_syn{index}_register.csv'.
    syn : Synapse instance, optional
        Object of the class Synapse, containing the required index in its attribute indexes['syn'].
        Necessary if the options 'syn' or 'resp' are used.
        Default value: None.
    **kwargs :
        Target values for filtering the register.
        The keys must belong to the columns of the table synapses_register:
        n_sites, area_tot_sites, res, border_intra, border_extra.
        See the attribute indentities['syn'] for details.

    Returns
    -------
    found : list of int
        References (indexes) of the candidate saves, which feature the desired properties.
        The list is empty if no save was found latching the criteria.
    """
    register = load_register(option, syn)
    if type(register) == pd.core.frame.DataFrame: # if a register exists
        filters = [key for key in kwargs.keys() if key in register.columns and key != 'Index']
        found = [i for i,row in register.iterrows() if all(row[col]==kwargs[col] or (np.abs(row[col]-kwargs[col])<= 10e-11) for col in filters)]
        # to solve the floating precision problem with read_csv: np.abs(row[col]-kwargs[col])<= 10e-11
        return found
    else:
        return []

def identity_summary(option, obj):
   """Generates dictionaries summarizing the main synapse/response attributes.
   The keys correspond to the columns of the registers of saved synapses/responses.
   The interest is to ease their access in the functions for saving and retrieving data.

    Parameters
    ----------
    option : {'syn', 'resp'}
        * 'syn': to get the summary of the synapse's main attributes,
        * 'resp': to get the summary of the stimulation pattern under which was generated the response currently stored in resp.
    obj : {Synapse, Pattern}
        Object of classs Synapse is necessary with option 'syn'.
        Object of classs Pattern is necessary with option 'resp'.

    Returns
    -------
    identities : dict
        Keys with option 'syn':
        * 'dim'
        * 'n_sites'
        * 'area_tot_sites'
        * 'res'
        * 'border_intra'
        * 'border_extra'
        * 'N_AMPA'
        * 'E_AMPA'
        * 'V0'
        * 'R'
        * 'C'.
        Keys with option 'resp':
        * 'mode'
        * 'tstep'
        * 'dur'
        * 'nit'
        * 'start'
        * 'n_stim'
        * 'f_stim'
        * 'fenv'
        * 'fmax'
        * 'fmin'.
    """
    if option == 'syn':
        identity = {key:obj.__dict__[key] for key in ['dim', 'n_sites', 'area_tot_sites', 
                                                    'res', 'border_intra', 'border_extra',
                                                    'N_AMPA', 'G_AMPA', 'E_AMPA', 'V0', 'R', 'C']}
    if option == 'resp':
        identity = {key:obj.__dict__[key] for key in obj.__dict__.keys() if (key != 'params' and key != 'spat')}
        for key in obj.params :
            identity[key] = obj.params[key]
    return identity

def attribute_index(option, syn):
   """Finds if the current synapse/response has already been saved in the registers. 
   Attributes the corresponding index.

    Parameters
    ----------
    option : {'syn', 'resp'}
        * 'syn': to attribute an index to the synapse, from the register of all synapses.
        * 'resp': to attribute an index to the response stored in the attribute "resp", from the register of the responses of the synapse.
    syn : Synapse instance
        Object of the class Synapse.

    Returns
    -------
    index : int
        References of the synapse/response instance in the registers of saves, also pointing to their respective directories.
    """
    option_sup = {'syn':'all_syn', 'resp':'syn'}[option] # to get the register one level above
    found = filter_register(option_sup, syn, **syn.identities[option])
    if len(found) == 0:
        index = None
    else:
        index = found[0]
    return index

def register_instance(option, syn):
   """Appends the current synapse/response in the registers. Sets an index. 

    Parameters
    ----------
    option : {'syn', 'resp'}
        * 'syn': to register the synapse, in the register of all synapses.
        * 'resp': to register the response stored in the attribute "resp", in the register of the responses of the synapse.
    syn : Synapse instance
        Object of the class Synapse.

    Returns
    -------
    index : int
        References of the synapse/response instance in the registers of saves, also pointing to their respective directories.
    """
    # Creating a dictionary summarizing the current instance:
    new_row = syn.identities[option]
    # Checking if it has not been appended yet:
    index = attribute_index(option, syn)
    if index == None:
        option_sup = {'syn':'all_syn', 'resp':'syn'}[option] # to get the register one level above
        create_register(option_sup, syn)
        register = load_register(option_sup, syn)
        index = len(register.index) # set index, starting at 0 if first registration
        with open(path_register(option_sup, syn), 'a') as file:
            pd.DataFrame(new_row, index=[index]).to_csv(file, header=False, float_format='%.16g')
    return index

def create_coords_ref(coords) :
    """Definition of reference names to save the responses for the points where glutamate concentrations are recorded.

    Parameters
    ----------
    coords : array-like, dtype = int
        Array of dimension (4, 2), containing the locations of the points where glutamate concentrations are recorded.
        See compute_coords_measures() for details.

    Returns
    ----------
    list of str
        List containing the reference names for the points.
        A reference name is a string that can be used to specify filnames when saving the results.
        Each name is of the shape : 'x-y', with x and y the coordinates of the point.
    """
    return ['{}-{}'.format(coords[i,0], coords[i,1]) for i in range(coords.shape[0])]



if __name__ == "__main__":
    # Printing the documentation:
    print("""This module is not designed to be executed in isolation in the shell. Please import UBC_synapse.
        Below is displayed the documentation about this module. Press 'q' to escape.""")
    print(__doc__)
    # Retrieving functions defined in this script:
    def is_local(object):
        return isinstance(object, FunctionType) and object.__module__ == __name__
    func_list = [value for name, value in getmembers(sys.modules[__name__], predicate=is_local)]
    for func in func_list :
        print(help(func))