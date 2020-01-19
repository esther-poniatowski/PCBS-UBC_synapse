#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class Synapse, for building models of synapses between mossy fibers and UBCs."""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"

print("synapse successfully imported")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
try :
    from UBC_synapse.synapse_morph import create_all_matrixes
    from UBC_synapse.stimulation_sim import compute_coords_measures, execute_c_code
    from UBC_synapse.visualize import print_matrix
    from UBC_synapse.manage_data import *
except ModuleNotFoundError :
    print("This module is not designed to be executed in isolation in the shell. Please import UBC_synapse.")


class Synapse:
    """Synapse model, defined by:
        * morphological properties: matrixes representing the UBC glomerule and its individual synaptic sites,
        * functional properties: diffusion of glutamate neurotransmitters in the cleft, dynamics of AMPA receptors, variations of the UBC's membrane potential.

    Attributes
    ----------
    dim : int
        Diameter of the glomerulus (entire synapse), in px.
    n_sites : int
        Number of individual glutamate release sites.
    area_tot_sites : float
        Cumulative area of glutamate release, in m².
    res : float
        Sptatial resolution of the model , in m (conversion factor). 
    border_intra : int
        Margin in px between the most external synaptic site and the synapse border.
    border_extra : int
        Margin in px between the synapse border and the limit of the working matrix.
    I : array_like, dtype = int
        Matrix representing the glomerular area (entire synapse).
        It contains values "1" within a disk representing the synapse, and values "0" elsewhere.
    O : array_like, dtype = int
        Matrix representing the extraglomerular area in the working matrix. O is the mirror of I, obtained by 1 - I.
    S : array_like, dtype = int
        Matrix representing the individual synaptic sites areas.
        The synaptic sites are represented by values "1".
        The centers of the synaptic sites remain spotted by values "2".
        The other pixels contain values "0".
    M_site : array-like, dtype = int
        Small patch (matrix) representing a single synaptic site area.
    N_AMPA : int
        Number of AMPA channels per pixel of individual synaptic site.
    G_AMPA : float
        Conductance of AMPA channel.
    E_AMPA : float
        Inversion potential of AMPA channel.
    V0 : float
        Resting potential of the UBC's membrane.
    R : float
        Resistance of the UBC's membrane.
        R = 1/G_leak, with G_leak = leak conductance.
    C : float
        Capacitance of the UBC's membrane.
    resp : dict
        Dictionary containing a response of the synapse under a particular stimulation. 
        It can be the result of the last simulation run, or a response retrived from the saved data for the current synapse instance.
        Keys:
        * 'coords': array of dimension (n, 2), containing the locations of the points where glutamate concentrations are recorded ([row, col]),
        * 'glus': list of Series, each one being the evolution of the glutamate concentrations at a location in 'coords',
        * 'AMPAtot': Series, evolution of the activation state of the AMPA channels (mean activation in the population of channels over the synapse),
        * 'V': Series, evolution of the UBC's membrane potential.
    identities : dict of dict
        Dictionaries of the main synapse/response attributes, corresponding to the columns of the registers of saved synapses/responses.
        The interest is to ease their access in the functions for saving and retrieving data.
        Keys:
        * 'syn': summary of the synapse's main attributes.
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
        * 'resp': summary of the stimulation pattern under which was generated the response currently stored in resp.
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
    indexes : dict of int
        References of the synapse/response instance in the registers of saves, also pointing to their respective directories.
        Keys:
            * 'syn': index of the synapse in synapses_register.csv,
            * 'resp': index of the current response stored in resp in responses_synx_regisster.csv.
    
    Main methods
    ------------
    instantiate (class method) - Instantiates several synapses withh various diameters.
    simulate_stimulation - Runs a simulation with the synapse, provided a stimulation pattern.
    save_response - Saves the last response computed for the synapse.
    visualize_synapse - Displays the synapse morphology.
    
    See also
    ----------
    synapse_morph - Module for building the matrixes I, O, S, M_site.
    stimulation_sim - Module for simulating the behavior of a synapse in response to a stimulation pattern, and interfacing with the c_code, 
    stimulation_pattern - Class for building a stimulation pattern.
    manage_data - Module for saving and retrieving the results.
    """

    # ===========
    # Constructor
    # ===========

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs : 
            Arguments corresponding to the desired values of the attributes.
            See the function create_all_matrixes from the module create_matrix for details about arguments and computations.
            Note: For "dim", "n_sites", "area_tot_sites", the value provided by the user might not be the final attribute of the synapse instance.
            The final value stored in the attribute is adjusted to match several constraints.
            dim : int
                Default value: None.
            n_sites : int
                Default value: 37.
            area_tot_sites : float
                Default value: 25e-12 m² (between 12 and 40 um², in vivo references given by Mungi).
            res : float
                Default value: 0.2e-6 m, corresponding to 1 px <=> 0.2 um. 
            border_intra : int
                Default value: 2 px.
            border_extra : int
                Default value: 2 px.
            N_AMPA : int
                Default value: 3.
            G_AMPA : float
                Default value: 20e-12 Ohm-1 (~10 pS)
            E_AMPA : float
                Default value: 0 mV
            V0 : float
                Default value: -75e-3 (-75 mV)
            R : float
                Default value: 652e6 (~100 MOhm, Van Drop 2014)
            C : float
                Default value: 17e-12 F (~10 pF, Van Drop 2014)
            verbose : bool, optional
                If True, the updated/suggested parameters and information will be printed.
                Default value: False.
        """
        self.res = kwargs.get('res', 0.2e-6)
        self.border_intra = kwargs.get('border_intra', 2)
        self.border_extra = kwargs.get('border_extra', 2)
        self.S, self.I, self.O, self.M_site, self.n_sites, self.area_tot_sites, self.dim = create_all_matrixes(n_sites=kwargs.get('n_sites', 37), 
                                                verbose=kwargs.get('verbose', False),
                                                dim=kwargs.get('dim', None), 
                                                area_tot_sites=kwargs.get('area_tot_sites', 25e-12),
                                                res=self.res, border_intra=self.border_intra, border_extra=self.border_extra)
        self.N_AMPA = kwargs.get('N_AMPA', 3)        
        self.E_AMPA = kwargs.get('E_AMPA', 0)         
        self.G_AMPA = kwargs.get('G_AMPA', 20e-12)
        self.V0 = kwargs.get('V0', -75e-3)         
        self.R = kwargs.get('R', 652e6 )          
        self.C = kwargs.get('C', 17e-12)
        self.resp = {'coords': None,
                    'glus': None,
                    'AMPAtot': None,
                    'V': None}
        self.identities = {'syn': identity_summary('syn', self), 'resp':None}
        self.indexes = {'syn': None, 'resp':None} # initializing to None to enable top directory creations
        self.indexes['syn'] = attribute_index('syn', self) # setting index if already saved

    # =============
    # Class methods
    # =============
   
    @classmethod
    def instantiate(cls, dim_list, **kwargs):
        """Class method. Instanciates several synapses with various dimensions.

        Parameters
        ----------
        dim_list : list of int
            Target diameters in px for the distinct glomerules.
        **kwargs :
            Additional non default arguments, target values for other attributes:
            n_sites, area_tot_sites, res, border_intra, border_extra.
            See the documentation of the class Synapse for the definition of the **kwargs arguments.

        Returns
        -------
        syn_dict : dictionary of Synapse instances.
            Keys: diameters of the synapses
            Values: Synapses objects instanciated.
        """
        instances = (cls(dim=dim, **kwargs) for dim in dim_list)
        syn_dict = {str(inst.dim):inst for inst in instances}
        print("Synapses instanciated: {}".format(syn_dict.keys()))
        return syn_dict

    @classmethod
    def search_synapses(cls, **kwargs):
        """Class method. Finds the saved synapses which feature specific properties.

        Parameters
        ----------
        **kwargs :
            Criteria for filtering the synapses. 
            The keys must belong to the columns of the table synapses_register. 
            See the attribute indentities['syn'] for details.

        Returns
        -------
        found : list of int
            References (indexes) of the candidate synapses, which feature the desired properties.
        """
        found = filter_register('all_syn', **kwargs)
        print("Indexes of found synapses:{}".format(found))
        register = load_register('all_syn').iloc[found]
        return found, register

    @classmethod
    def retrieve_synapse(cls, synindex):
        """Class method. Loads a particular synapse which has been previousely saved.

        Parameters
        ----------
        synindex : int
            Reference of the synapse in the table synapses_register.

        Returns
        -------
        syn : Synapse object
        """
        try :
            synapses_register = load_register('all_syn')
            # Retrieving the desired attributed
            attrs = synapses_register[synapses_register["Index"]==synindex].drop(['Index'], axis=1).to_dict('records')[0]
            # Intanciation of the synapse:
            syn = cls(**attrs)
        except :
            print("No synapse saved under this index.")
            syn = None
        finally :
            return syn

    # ================
    # Instance methods
    # ================

    def search_responses(self, **kwargs):
        """Finds the saved responses which feature specific properties.

        Parameters
        ----------
        **kwargs :
            Criteria for filtering the responses. 
            The keys must belong to the columns of the table responses_synx_register.
            See the attribute indentities['resp'] for details.
        
        Returns
        -------
        found : list of int
            References (indexes) of the candidate responses, which feature desired properties.
        """
        found = filter_register('syn', self, **kwargs)
        print("Indexes of found responses for this synapse:{}".format(found))
        register = load_register('syn', self).iloc[found]
        return found, register

    def retrieve_response(self, respindex):
        """Loads a particular saved response.
        It fills the dictionary attribute "resp" with different components of the simulation.
        See the attribute "resp" for more details.
        
        Parameters
        ----------
        respindex : int
            Reference of the response in the table responses_synx_register.
        """
        try :
            # Unpacking the data in the .resp attribute:
            self.indexes['resp'] = respindex
            path_dir = path_directory('resp', self)
            self.resp['coords'] = pd.read_csv(os.path.join(path_dir, 'coords.csv'), index_col=0).to_numpy()
            coords_ref = create_coords_ref(self.resp['coords'])
            self.resp['glus'] = tuple(pd.read_csv(os.path.join(path_dir, 'resglu{}.csv'.format(ref)), index_col=0) for ref in coords_ref)
            self.resp['AMPAtot'] = pd.read_csv(os.path.join(path_dir, 'resAMPAtot.csv'), index_col=0)
            self.resp['V'] = pd.read_csv(os.path.join(path_dir, 'resV.csv'), index_col=0)
            print("Response retrieved.")
        except FileNotFoundError:
            print("No response saved under this index.")

    def simulate_stimulation(self, patt):
        """Runs a simulation with the current instance.
        If the simulation has already been run previousely, the response is retrived in the attribute "resp".

        Parameters
        ----------
        patt : Pattern object
            Pattern of stimulation, under which the response will be recorded.
        """
        # Defining the response:
        self.identities['resp'] = identity_summary('resp', patt)
        respindex = attribute_index('resp', self)
        # Running the simulation if no response has been computed for this pattern:
        if respindex == None :
            print('Running the simulation. It may take some time.')
            self.resp['coords'], self.resp['glus'], self.resp['AMPAtot'], self.resp['V'] = execute_c_code(self, patt)
            print("Simulation completed.")
        # Retrieving the existing response otherwise:
        else:
            print("Response already computed.")
            self.retrieve_response(respindex)

    def save_response(self):
        """Saves the response of the current instance, stored in the attribute "resp".
        It registers the synapse in the register of all synapses, if it has never been saved before.
        It registers the response in the register of responses of the synapse.
        It creates a sub-folder for the response in the directory of the synapse, to drop the data in .csv files.
        """
        self.indexes['resp'] = attribute_index('resp', self)
        # Checking if the attribute "resp" is not empty:
        if not type(self.resp['coords']) == np.ndarray:
            print("Response is empty. Please run a simulation.")
        # Checking if the target response has already been registered:
        elif self.indexes['resp'] == None:
            # Registering the synapse if necessary:
            self.indexes['syn'] = register_instance('syn', self)
            # Registering the response and setting its new index:
            self.indexes['resp'] = register_instance('resp', self)
            create_directory('resp', self)
            # Exporting the contents of the attribute "resp" to csv files:
            path_dir = path_directory('resp', self)
            coords_ref = create_coords_ref(self.resp['coords'])
            pd.DataFrame(self.resp['coords']).to_csv(os.path.join(path_dir, 'coords.csv'))
            for i in range(len(coords_ref)) :
                self.resp['glus'][i].to_csv(os.path.join(path_dir, 'resglu{}.csv'.format(coords_ref[i])), header=True)
            self.resp['AMPAtot'].to_csv(os.path.join(path_dir, 'resAMPAtot.csv'), header=True)
            self.resp['V'].to_csv(os.path.join(path_dir, 'resV.csv'), header=True)
            print("Saved: response at index {} for synapse {}.".format(self.indexes['resp'], self.indexes['syn']))
        else:
            print("Response already registered at index {} for synapse {}.".format(self.indexes['resp'], self.indexes['syn']))


    def visualize_synapse(self, print_coords=False):
        """Displays the matrixes representing the synapse.
 
        Parameters
        ----------
        print_coords : bool, optional
            If True, the points where glutamate concentrations were recorded will be displayed.

        See also
        --------
        print_matrix in the module visualize
        """
        S = self.S.copy()
        I = self.I.copy()
        if print_coords :
            print_matrix(S, I, self.res, print_coords=True, coords=self.resp['coords'])
        else :
            print_matrix(S, I, self.res)




if __name__ == "__main__":
    # Printing the documentation:
    print("Below is displayed the documentation about this module. Press 'q' to escape.")
    print(__doc__)
    print(help(Synapse))