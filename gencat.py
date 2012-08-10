import numpy as np
from math import *
import scipy as sp
import scipy.interpolate
from parameters import *

#Functions.............................
def err_magnitude():

	#Computing photons number
	N_sig = (3631.0 * 1.51 * pow(10, 7) * pow(10,-0.4 * m) * in_r[i] * t_exp[i] * tel_surface) / n_pix #object photons per pixel
	N_sky = in_sky[i] * t_exp[i] * tel_surface * pix_size #sky photons per pixel

	#Computing Signal to Noise (StoN)
	StoN = (sqrt(n_pix * n_exp) * N_sig) / sqrt(N_sig + N_sky + RN * RN)
	noise_ctn = 2.5 * log10(1 + 0.02)
	err_m_obs = 2.5 * log10(1 + 1 / StoN)
	err_m_obs = sqrt((err_m_obs * err_m_obs) + (noise_ctn * noise_ctn))

	return err_m_obs 

#Loading input catalog file.............
cat = np.loadtxt(incat_file + incat_file_extention, usecols = (index_id, index_m, index_z, index_t), unpack = 'True')

id = cat[0]
m_ref = cat[1]
z = cat[2]
t = cat[3]
n_gal = len(m_ref)

#Loading exposure times................
t_exp = np.loadtxt(exp_t_file, unpack = True)

#Loading filters....................... 
filt_list = np.loadtxt(filt_folder + filt_list, dtype = 'string', unpack = True)
n_filt = len(filt_list)
filt_list = np.hstack(([ref_filt], filt_list))

#Computing in_r.........................
y_r = {}
x = {}
in_r = np.zeros(n_filt + 1)
for i in range(n_filt + 1): 
	filt = filt_list[i]
	R = np.loadtxt(filt_folder + filt, unpack = True)

	#Defining lambda range
	x[filt] = np.arange(R[0].min(), R[0].max(), dx)

	#Interpolation
	y_r[filt] = sp.interpolate.interp1d(R[0], R[1])(x[filt]) 
	in_r[i] = (y_r[filt] / x[filt]).sum() * dx

in_r0 = in_r[0]
in_r = in_r[1:]

#Computing in_sky.......................
SKY = np.loadtxt(sky_file, unpack = True) 
in_sky = np.zeros(n_filt + 1)
for i in range(n_filt + 1):
	y_sky = sp.interpolate.interp1d(SKY[0], SKY[1])(x[filt])
	in_sky[i] = (y_sky * y_r[filt]).sum() * dx

in_sky = in_sky[1:]

#Loading seds...........................
sed_list = np.loadtxt(sed_folder + sed_list, dtype = 'string', unpack = True)
n_sed = len(sed_list) 

#dt = (t.max() - t.min()) / (n_sed - 1)
#t_range = np.arange(t.min(), t.max() + dt, dt )
dt = 1
t_range = np.arange(0, 66)

s = {}
for sed in sed_list: s[sed] = np.loadtxt(sed_folder + sed, unpack = True)

f_noiseless = open(f_noiseless_file, "w") 
f_noisely = open(f_noisely_file, "w") 

#Iterating over all the galaxies
for n in range(n_gal):
	f_noiseless.write("%d " % id[n])
	f_noisely.write("%d " % id[n])

	s_ind_high = np.searchsorted(t_range, t[n])
	s_high = sed_list[s_ind_high]
	s_low = sed_list[s_ind_high - 1]
	coef_high = (t_range[s_ind_high] - t[n]) / dt 
	coef_low = 1 - coef_high  

	#Computing in_s................
	in_s = np.zeros(n_filt + 1)
	for i in range(n_filt + 1):
		filt = filt_list[i]
		y_s_low = sp.interpolate.interp1d(s[s_low][0] * (1 + z[n]), s[s_low][1])(x[filt]) 
		y_s_high = sp.interpolate.interp1d(s[s_low][0] * (1 + z[n]), s[s_low][1])(x[filt]) 
		y_s = coef_high * y_s_high + coef_low * y_s_low
		in_s[i] = (y_s * y_r[filt] * x[filt]).sum() * dx
		
	in_s0 = in_s[0]
	in_s = in_s[1:]

	#Generating magnitudes.........
	for i in range(n_filt):
		m = m_ref[n] + 2.5 * (log10(in_s0) + log10(in_r[i]) - log10(in_r0) - log10(in_s[i])) 
				
		err_m_obs = err_magnitude()
		rn = np.random.normal()
		m_obs = m + rn * err_m_obs	

		f_noiseless.write("%4.4f " % m)
		f_noisely.write("%4.4f %4.4f " % (m_obs, err_m_obs))

	f_noiseless.write("%4.4f %2.2f\n" % (z[n], t[n]))
	f_noisely.write("%4.4f %2.2f\n" % (z[n], t[n]))

f_noiseless.close()
f_noisely.close()
