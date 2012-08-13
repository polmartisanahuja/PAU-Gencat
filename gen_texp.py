import numpy as np
from parameters import *

tray_matrix = np.loadtxt(filt_folder + tray_matrix_file, dtype = 'string')
filt_names = np.loadtxt(filt_folder + filt_names_file, dtype = 'string')

n_tray = len(texp_tray)
n_filt = len(filt_names)

#Creating the exposure times matrix.....................
texp_matrix = np.ones((6,16))

for i in range(n_tray):
	texp_matrix[i][:8] = texp_tray[i] * 1. #Central CCDs
	texp_matrix[i][8:14] = texp_tray[i] * 0.75 #Peripheral CCDs
	texp_matrix[i][14:16] = texp_tray[i] * 0.375 #Corner CCDs

#Creating the exposure time list......................
texp = np.zeros(n_filt)
filt_list = {}
for filt in filt_names: filt_list[filt] = 0 
for i in range(n_tray):
	for j in range(16):
		filt_list[tray_matrix[i][j]] += texp_matrix[i][j]
for i in range(n_filt): texp[i] = filt_list[filt_names[i]]
#np.savetxt(texp_file, texp, fmt = "%6.6f")
