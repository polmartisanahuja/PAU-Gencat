import numpy as np
from pau_parameters_gencat import *

filt_folder, filt_names = np.loadtxt(filt_path + filt_names_file, dtype = 'string', unpack = True)

n_tray = len(texp_tray)
n_filt = len(filt_names)
n_BB = len(texp_BB)
print "n_tray=", n_tray
print "n_filt=", n_filt
print "n_BB=", n_BB 

#Creating the exposure times matrix.....................
texp = np.zeros(n_filt)

for i in range(n_BB): texp[i] = texp_BB[filt_names[i]] 
for i in range(n_filt - n_BB): texp[n_BB + i] = texp_tray[i/n_filt_tray]

#Creating the exposure time list......................
np.savetxt(filt_path + texp_file, texp, fmt = "%6.6f")
