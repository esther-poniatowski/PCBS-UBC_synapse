# PCBS-UBC_synapse


### Modeling signal filtering in the synapse between Mossy Fibers and Unipolar Brush Cells
---

The cerebellum is a brain area involved in the representation of movement and position. This high-level representation seems to be achieved through the integration and the disambiguation of sensory signals. This idea motivates to investigate the determinants of the signal processing, that might lie at the level of the neurons.  

In the cerebellum, the signal flows through a feedforward circuit :  
(1) It receives direct sensory inputs from the vestibular system of the inner ears, conveyed by the *Mossy Fibers (MFs)*. This signal is produced by sensors of inertial translational and rotational acceleration.  
(2) In the cerebellum, the mossy fibers contact the *Unipolar Brush Cells (UBCs)*.  
![Feedforward network in the cerebellum](https://github.com/esther-poniatowski/PCBS-UBC_synapse/blob/master/UBC_references/fig1_cerebelllum_network.PNG)

The cerebellum features a high density of these remarkable interneurons, which display a very particular type of synapses. Indeed, the many convolutions of the UBCs' synaptic clefts lead to a giant contact surfaces, riddled with neurotransmitter release sites.  
This observation triggers the hypothesis that the synapse MF-UBC could potentially play a role in the transformation of the sensory signal encoding acceleration. More precisely, the synapse could act like a filter, leading to the temporal integration of the signal. This process would be a step to convert the sensory detection of acceleration into the cognitive perception of movement.

**The goal of the project is to explore the influence of the synapse's morphological and physiological properties on the signal transformation.**


### Usage
---

Please use the Jupyter notebook provided. This demonstration is designed to guide the user in the usage of the project, from building of the model to visualisating of the results.

The dependencies required for running the project are standard libraries from Python 3. They are nonetheless listed in several files in the folder `UBC_installation/` if necessary. They could be installed either via conda or virtualenv, thanks to the files environment.yml or requirements.txt.


### Progammation objectives
---

This project was drawn from an original project conducted in 2018-2019 in the context of an biology L3 intership, under the direction of Guillaume Dugué. The focus was directed, as priority, on investigating hypotheses about the synapse functions. Therefore, the codes I produced were hardly readable, all the more since carrying on this project was an opportunity for me to learn and practice programmation.
For the current project, aim rather consisted in providing a more sharable and automated tool, in view of further investigation based on the model. To get in shape the project, the main work leads were the following :  
* Translating the code from R language to Python language.
* Improving the code and the structure of the project, in terms of efficiency and understandability.
* Presenting the results in a notebook for simplifying the use.

I tried to implement was I had learned in the context of PCBS and from other programming resources:
* **Documentation** - I leaned about NumpyDocstyle syntax. I provided as much details as possible in all my modules, classes, and functions.
* **Don't Repeat Yourself** principle and Separation of Concerns. - I tried to make the code more flexible and to automate each step (search of parameters, instantiation of models, running simulations, saving and loading data, visualizing results). For this purpose, I rethought the architecture of the code, focusing on modularity: implementing the steps in different modules, dividing the unitary tasks in individual functions. It led me to start back from scratch on many aspects of the project. Namely, I designed a new method to build the matrixes representing the synapses, in order to optimize automatically the geometric constraints. I also worked on the structure of the files' tree to store data more efficiently. I tried to find conventions to keep track and refer to saved data, so that it could be browsed and accessed more easily.
* **Oriented Object Programming** - Forthe first time, I wanted to put in practice Oriented Object Programming. In this view, I thought that a modelling project was more adapted than data analysis. I built two classes to represent the synapses and the stimulation patterns. Those two classes consitute the functional interface for the user: they gather the main funcitonalities, without requiring to dig into the code details of other modules.
* **Packaging** - I readed about virtual environments and packages distribution. Even if my project is not meant to be distributed via Anaconda Cloud or PiPY, I tried to get the project in shape so that it could be shared. I built a "cookiecutter" template for my future projects. I listed the dependencies of my project in setup files.


### Statement of accounts
---

The main challenges I met were the following:
* **Compatibility between C and Python** - Optimizing running times required that the computations were carried out by a C code. It implied to call the C code from the Python code, and to retrieve the outputs yielded by the C code. To do so, the object types defined in Python have to fit specific formats, secified thanks to the library ctypes. Fixing the compatibility problems was time consuming. However, it was the opportunity for to deepen my knowledge about this technical aspect of programming.
* **Organizing the project structure** - I hesitated before settling on a stable structure. The most challenging was to find a procedure for saving the data generated by the codes. I opted for a separation in different folders according to the synapses, and a further division in sub-folders according to the simulations. I chose to keep track of the saved data in register files (one for all the synapses, other ones for the simulations for each synapse). I built functions to scan the registers, and recover the objets if needed. Nevertheless, this approach turned out to complexify the codes.

Lines of improvement:
* **Sphinx** - Owing to the multiplication of modules and functions, the user might be confused in the usage of the project, in spite of the Jupyter notebook. It might be useful to gather the documentation in a clearer shape, using Sphinx. However, I had not enough time to dedicate to learning how to build such a documentation.
* **Tests** - Initially, I had planned to implement a test-driven approach of writing code, using pytest. Nevetheless, I did not systematically have the reflex to write case tests in parallel with building code. I kept the tendency to carry out unit tests in a Jupyter notebook
* **Data generation** - In spite of using a C code, the running times remain long. For optimizing times, I had planned to implement parallel computing, with the module Dask. Since I lacked time to put it in practice in the context of the project, the examples provided in the Jupyter notebook are scarce and limited. Their purpose is to illutrate briefly the usage of the project, rather than analysing the model in depth, which would require multiple simulations on long time series. I also wanted to write functions for downsampling the time series, and performing analyses (Fourier, derivatives...).
* **Visualization of the results** - I would like to make more modular functions for visualizing each component of the responses. They could be called individually to build a facetted plot.

To conclude, I yearn to pursue this project by appending new functionalities. Having taken advantage of PCBS for building a more structured and documented project may facilitate further development of the model.


### References
---

Tsayem Ngueguim Idriss, Mémoire, *Bases cellulaires du traitement de l’information vestibulaire dans le cervelet: mise en place d’un modèle de fonctionnement synaptique et d’une nouvelle approche pour l’observation anatomique de l’origine des fibres vestibulaires aférentes*, 2018.

I. M. Raman, L. O. Trussell, *The mechanism of a-amino-3-hydroxy-5-methyl-4-isoxazolepropionate receptor desensitization after removal of glutamate*, Biophys. J.,68(1):137–146, 1995.

S. van Dorp, C. I. De Zeeuw, *Variable timing of synaptic transmission in cerebellar unipolar brush cells*, Proceedings of the National Academy of Sciences, 111(14):5403–5408, 2014.

S. van Dorp, C. I. De Zeeuw, *Forward Signaling by Unipolar Brush Cells in the Mouse Cerebellum*, Cerebellum, 2015.

V. Zampini, J. K Liu, M. A. Diana, P. P. Maldonado, N. Brunel, S. Dieudonne, *Mechanisms and functional roles of glutamatergic synapse diversity in a cerebellar circuit*, eLIFE, 2010.

E. Mugnainia, G. Sekerkováa, M. Martina, *The unipolar brush cell: A remarkable neuron finally receiving deserved attention*, Brain Research Reviews, 2010.

T. S. Balmer, L. O. Trussell, *Selective targeting of unipolar brush cell subtypes by cerebellar mossy fibers*, bioRxiv, 2019.

M. Masugi-Tokita, E. Tarusawa, M. Watanabe, E. Molnar, K. Fujimoto, R. Shigemoto, *Number and Density of AMPA Receptors in Individual Synapses in the Rat Cerebellum as Revealed by SDS-Digested Freeze-Fracture Replica Labeling*, The Journal of Neuroscience, 2007.

T. Caitlin Smith, L. Y. Wang,2 J. R. Howe, *Heterogeneous Conductance Levels of Native AMPA Receptors*, The Journal of Neuroscience, 2000.

