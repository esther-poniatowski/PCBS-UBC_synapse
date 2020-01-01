#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Modeling of synapses between mossy fibers and unipolar brush cells (UBCs) in the cerebellum.

The goal of this module is to investigate the properties of a specific type of synapses in the cerebellum: between mossy fibers and unipolar brush cells (UBCs).
These synapses are interesting because they feature gigantic dimensions, with numerous neurotransmitter (glutamate) release sites.
They could potentially play a role in the signal transmission in the cerebellum, acting like filters. They could be involved in the temporal integration of the signal, from the sensory detection of acceleration, to the cognitive perception of movement. 

The aim is to assess the effect of the diameter of the synapse on the signal transmission.
For this purpose, the package achieves the following tasks:
    1. Modeling the morphology of synapses.
    2. Simulating their behavior in response to a stimulation pattern (events of glutamate release).

Usage
-------
The Jupyter notebook provided in the package was designed to guide the user in the exploration of the functionalities.

Structure of the package
------------------------
PCBS_UBC_synapse/
│
├── UBC_simulation.ipynb    <- notebook to guide the user into the project
├── README.md               
├── LICENSE.txt
├── AUTHORS.rst
├── environment.yml         <- list of dependencies, file usable by conda
├── requirements.txt        <- list of dependencies, file usable by virtualenv
├── setup.py                <- configuration file for an enventual distribution via PiPy
├── MANIFEST.in             <- specification of the non-.py files to be included in the project
├── UBC_synapse/            <- source codes directory
│   └── __init__.py
│   └── create_matrix.py
│   └── synapse.py
│   └── patterns_release.py
│   └── simulate.py
│   └── visuaize.py
│   └── tools.py
│   └── c_code/
│       └── run_simulation.c
│       └── run_simulation.so
├── UBC_data/               <- directory for storing the data generated by running the codes
│   └── register_syn.csv
│   └── syn1/
│       └── syn1.txt
│       └── register_responses.csv
│       └── patt1/
│           └── resgluxx-xx.csv
│           └── resAMPAtot.csv
│           └── resV.csv

Modules
-------
All the Python modules are located in the directory `UBC_synapse/`. The C library is contained in the sub-folder `c_codes/`.
    * create_matrix - Functions for building the morphology of a synapse.
    * synapse - Class for building synapses objects. 
    * patterns_release - Class for generating stimulation patterns, under which the synapses' behaviors will be observed.
    * run_simulation (C code) - Core computations for simulating the response of a synapse under a stimulation pattern.
    * simulate - Functions for interfacing with the C code. 
    * visuaize - Functions for generating representations of the synapses' morphology and behaviors.
    * tools - General tools for other modules' computations.
A more detailed documentation for individual modules is provided in each one, with their respective functions.

Data
-------
The synapses and the data generated by running the package can be saved in the directory `UBC_data/`, thanks to the methods of the class Synapse.
The file `register_syn.csv` keeps track of the saved synapses: 
    * the different instances are classified according to their attributes,
    * each saved synapse is referenced by an index (number).
Each saved synapse is associated to a sub-directory named `sny{index}`. One directory contains:
    * a file `sny{index}` created with the library pickle, allowing to retrieve the object with its attributes,
    * a file `register_responses.csv`, classifying the stimulation trial with this synapse, and referencing them by an index (number),
    * for each stimulation trial, a sub-folder `patt{index}/` containing the time series of the response:
        * `coords.csv`: locations of the points on the synapse where the glutamate concentrations were recorded,
        * `resgluxx-xx.csv`: glutamate concentrations at location xx-xx (row-column),
        * `resAMPAtot.csv`: activation state of the post-synaptic AMPA receptors,
        * `resV.csv`: UBC membrane potential.
"""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"

print("UBC_synapse importation:")