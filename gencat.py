import numpy as np
from math import *
import scipy as sp
import scipy.interpolate
from parameters import *

#Loading input catalog file.............
cat = np.loadtxt(incat_file + incat_file_extention, usecols = (index_id, index_m, index_z, index_t), unpack = 'True')

id = cat[0]
m_ref = cat[1]
z = cat[2]
t = cat[3]
n_gal = len(m_ref)

#Loading filters & computing in_r......
filt_list = np.loadtxt(filt_folder + filt_list, dtype = 'string', unpack = True)
filt_list = np.hstack(([ref_filt], filt_list))
n_filt = len(filt_list)

y_r = {}
x = {}
in_r = np.zeros(n_filt)
for i in range(n_filt): 
	filt = filt_list[i]
	R = np.loadtxt(filt_folder + filt, unpack = True)

	#Defining lambda range
	x[filt] = np.arange(R[0].min(), R[0].max(), dx)

	#Interpolation
	y_r[filt] = sp.interpolate.interp1d(R[0], R[1])(x[filt]) 
	in_r[i] = (y_r[filt] / x[filt]).sum() * dx

#Loading templates.....................
sed_list = np.loadtxt(sed_folder + sed_list, dtype = 'string', unpack = True)
n_sed = len(sed_list) 
dt = (t.max() - t.min()) / (n_sed - 1)
t_range = np.arange(t.min(), t.max() + dt, dt )

s = {}
for sed in sed_list: s[sed] = np.loadtxt(sed_folder + sed, unpack = True)

#Generating magnitudes.................
f_noiseless = open(f_noiseless_file, "w") 
for n in range(n_gal):
	f_noiseless.write("%d " % id[n])

	s_ind_high = np.searchsorted(t_range, t[n])
	s_high = sed_list[s_ind_high]
	s_low = sed_list[s_ind_high - 1]
	coef_high = (t_range[s_ind_high] - t[n]) / dt 
	coef_low = 1 - coef_high  

	in_s = np.zeros(n_filt)
	for i in range(n_filt):
		filt = filt_list[i]
		y_s_low = sp.interpolate.interp1d(s[s_low][0] * (1 + z[n]), s[s_low][1])(x[filt]) 
		y_s_high = sp.interpolate.interp1d(s[s_low][0] * (1 + z[n]), s[s_low][1])(x[filt]) 
		y_s = coef_high * y_s_high + coef_low * y_s_low
		in_s[i] = (y_s * y_r[filt] * x[filt]).sum() * dx
		
		m = m_ref[n] + 2.5 * (log10(in_s[0]) + log10(in_r[i]) - log10(in_r[0]) - log10(in_s[i])) 

		if( (i != 0) & (i != (n_filt - 1))): f_noiseless.write("%4.4f " % m)
	f_noiseless.write("%4.4f %2.2f\n" % (z[n], t[n]))

f_noiseless.close()
