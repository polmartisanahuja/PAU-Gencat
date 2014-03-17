import sys
import numpy as np
from math import *
import scipy as sp
import scipy.interpolate
from pau_parameters_gencat import *

#Loading filters....................... 
filt_folder, filt_list = np.loadtxt(filt_path + filt_names_file, dtype = 'string', unpack = True)
n_filt = len(filt_list)

#Computing in_r.........................
y_r = {}
x = {}
in_r = np.zeros(n_filt)
for i in range(n_filt): 
	R = np.loadtxt(filt_path + filt_folder[i] + filt_list[i] + ".res", unpack = True)

	filt = filt_list[i]

	#Defining lambda range
	x[filt] = np.arange(R[0].min(), R[0].max(), dx)

	#Interpolation
	y_r[filt] = sp.interpolate.interp1d(R[0], R[1])(x[filt]) 
	in_r[i] = (y_r[filt] / x[filt]).sum() * dx

#Computing mean Lambda...................
mean_lam = np.zeros(n_filt)
for i in range(n_filt): 
	lam, r = np.loadtxt(filt_path + filt_folder[i] + filt_list[i] + ".res", unpack = True)
	mean_lam[i] = (lam * r).sum() / r.sum()	

#Loading exposure times................
texp = np.loadtxt(texp_file, unpack = True)
texp *= n_exp 

#Computing in_sky......................
SKY = np.loadtxt(sky_file, unpack = True) 
in_sky = np.zeros(n_filt)
for i in range(n_filt):
	filt = filt_list[i]
	y_sky = sp.interpolate.interp1d(SKY[0], SKY[1])(x[filt])
	in_sky[i] = (y_sky * y_r[filt]).sum() * dx
N_sky = in_sky * texp * tel_surface * pix_size #sky photons per pixel

#Computing m_lim.........................
A = (3631.0 * 1.51 * pow(10, 7) * in_r * texp * tel_surface) / n_pix
B = 5. / np.sqrt(n_pix)
C = RN * RN + N_sky

m_lim = -2.5 * np.log10((B / (2 * A)) * (B + np.sqrt(B*B + 4 * C)))

#np.savetxt(filt_path + 'PAU_lim_mag.txt', np.array([mean_lam,m_lim]).T, fmt = "%6.6f")
np.savetxt(filt_path + 'PAU_x2_lim_mag.txt', np.array([mean_lam,m_lim]).T, fmt = "%6.6f")
