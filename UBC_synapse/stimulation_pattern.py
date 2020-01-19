#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for generating stimulation patterns, under which the synapses' behaviors will be observed."""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"

print("patterns_release successfully imported")

import numpy as np
try :
    from UBC_synapse.visualize import print_pattern
except ModuleNotFoundError :
    print("""This module is not designed to be executed in isolation in the shell. Please import UBC_synapse.
            Below is displayed the documentation about this module. Press 'q' to escape.""")


def convert_time_to_steps(dur, tstep=1e-7):
    """Conversion of a duration in seconds to a number of time steps in the model.

    Parameters
    ----------
    dur : float
        Duration in seconds to be converted.
    tstep : float, optional
        Temporal resolution of the model, in s. Conversion factor, specifying the duration of 1 time step (iteration) in the model.
        Default: 1e-7 s, corresponding to 1 time step <=> 1e-7 s = 0.1 m.
    
    Returns
    -------
    int
        Number of time steps (iterations) corresponding to the given duration in the model (inferior integrer part).
    """
    return int(dur/tstep)


class Pattern:
    """Sitmulation pattern under which the behavior of the synapses will be observed.

    A pattern describes, in time, the events of glutamate release by the presynaptic neuron.

    Attributes
    ----------
    tstep : float
        Temporal resolution used for the pattern instance, in s. Conversion factor, specifying the duration of 1 time step (iteration) in the model.
    dur : float
        Duration of the pattern, in seconds.
    nit : int
        Number of iterations to run through the pattern, defining the length of the series.
    mode : {'stim', 'oscill'}
        Type of pattern.
            * 'stim': several release events at the beginning of the pattern, at a fixed frequency ;
            * 'oscill': continuous firing, whose frequency follows a sinusoidal variation.
    spat : array_like, dtype = int
        1D boolean array, encoding the release events by values '1'.
    params : dictionary
        Other parameters to define the pattern.
        See pattern_stim and pattern_oscill for details.
    """
    
    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs : 
            Arguments corresponding to the attributes.
            tstep : int, optional
                Defaultvalue: 1e-7 s, corresponding to 1 time step <=> 1e-7 s = 0.1 m.
            dur : float, optional
            nit : int, optional
                Default value: 20000.
            mode : {'stim', 'oscill'}, optional
                Default value: 'stim'.
        The other parameters passed by the user are stored in the argument params.
        See pattern_stim and pattern_oscill for details.
        """
        self.mode = kwargs.get('mode', 'stim')
        self.tstep = kwargs.get('tstep', 1e-7)
        if self.mode=='stim':
            self.dur = kwargs.get('dur', kwargs.get('nit', 20000)*self.tstep) # if not provided, conversion of nit (provided or default)
            self.nit = kwargs.get('nit', convert_time_to_steps(self.dur, self.tstep)) # conversion of dur, to ensure consistency
            self.params = {'start':kwargs.get('start', 1000), 
                            'n_stim':kwargs.get('n_stim', 1), 
                            'f_stim':kwargs.get('f_stim', 100)}
            self.pattern_stim()
        if self.mode=='oscill':
            self.dur = kwargs.get('dur', 1)
            self.nit = kwargs.get('nit', convert_time_to_steps(self.dur, self.tstep))
            self.params = {'fenv':kwargs.get('fenv', 10), 
                            'fmin':kwargs.get('fmin', 40), 
                            'fmax':kwargs.get('fmax', 150)}
            self.pattern_oscill()

    def pattern_stim(self):
        """Creates a stimulation pattern of type 'stim', in order to set the attribute spat.

        Parameters
        ----------
        The parameters are implicitly stored in the argument params:
            start : int
                Number of iterations before the first glutamate release event.
                Default value: 1000.
            n_stim : int
                Number of release events.
                Default value: 1
            f_stim : float
                Frequency of release, in Hz, determining the time interval between each event.
                Default value: 100 Hz.
            """
        start = self.params['start']
        n_stim = self.params['n_stim']
        f_stim = self.params['f_stim']
        interv = convert_time_to_steps(1/f_stim, self.tstep)
        spat = np.zeros(self.nit)
        for i in np.arange(start, n_stim*interv+1, step=interv):
            spat[i] = 1
        self.spat = spat.astype(np.bool_)

    def pattern_oscill(self):
        """Creates a stimulation pattern of type 'oscill', in order to set the attribute spat.

        Parameters
        ----------
        The parameters are implicitly stored in the argument params:
            fenv : float
                Frequency of global oscillations, i.e. of the "enveloppe" of the pattern.
                Default value: 10 Hz.
            fmin : float
                Minimum local frequency of firing, i.e. between two spikes at one given moment in the sequence.
                Default value: 40 Hz
            fmax : float
                Maximum local frequency of firing, i.e. between two spikes at one given moment in the sequence..
                Default value: 150 Hz.
        """
        fenv = self.params['fenv']
        fmin = self.params['fmin']
        fmax = self.params['fmax']
        a = 0.5*(fmax - fmin)
        b = 0.5*(fmax + fmin)
        i = 0
        spat = np.zeros(self.nit)
        while i < self.nit :
            time_interv = 1/(b+a*np.sin(2*np.pi*fenv*i*self.tstep)) # new delay before next spike, in s
            i = i + np.int(time_interv/self.tstep) # index of next spike 
            if i < self.nit :
                spat[i]=1
        self.spat = spat.astype(np.bool_)

    def visualize_patt(self):
        """Method for displaying the stimulation pattern.
        See print_patterin visualize for details."""
        print_pattern(self.spat, self.tstep)




if __name__ == "__main__":
    # Printing the documentation:
    print(__doc__)
    print(help(Pattern))