#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Parameters of the model."""

__authors__ = ("esther-poniatowski")
__contact__ = ("esther.poniatowski@ens.fr")
__version__ = "1.0.0"
__date__ = "01/01/2020"

# ==================
# Spatial parametres
# ==================
res = 0.2e-6	     # resolution : 1 px <=> 0.2 um
area_syn_tot = 40e-12 # Mungi
border_intra = 2     # margin in px between the most external sites and the synapse border
border_extra = 2     # margin in px between the synapse border and the limit of the working matrix

# ===================
# Temporal parameters
# ===================
tstep = 1e-7        # 1 time step <=> 1e-7 s = 0.1 ms
dur = 1             # pattern default duration

# =====================
# Functional parameters
# =====================
quant = 4           # quantum of glutamate release (mM)  : 4 mmol = 4000 umol
Diff = 7.6e-10      # diffusion coefficient (m2/s),
qin = Diff*tstep/(res**2) # equivalent of diffusion coefficient, adapted in px and tstep for the model 
qout = qin          # no diffusion barrier at the limit of glomerulus
N_AMPA_tot = 30     # number of AMPA channels per synapse
N_AMPA = 3          # number of AMPA channels per px :(10 px per synapse => 3 channels par px)
E_AMPA = 0          # inversion potential of AMPA channel (0 mV)
G_AMPA = 20e-12     # conductance of AMPA channel (~ 10 pS)
V0 = -75e-3         # resting potential of the membrane (-75 mV)
R = 652e6           # resistance of the membrane (~ 100 MOhm, Van Drop 2014), R = 1/G_l = leak conductance
C = 17e-12          # capacitance ofthe membrane (~ 10 pF, Van Drop 2014)

# ================================
# Dictionary of default parameters
# ================================
param_dict = {
    'dim':82, 'n_sites':55, 'n_px_site':10,
    'dim_ref':82, 'border_intra':3, 'border_extra':3, 'res':res,
    'quant':quant, 'qin':qin, 'qout':qout, 
    'N_AMPA':N_AMPA, 'E_AMPA':E_AMPA, 'G_AMPA':G_AMPA, 
    'R':R, 'C':C, 'V0':V0}


if __name__ == '__main__':
	print("parameters successfully imported")