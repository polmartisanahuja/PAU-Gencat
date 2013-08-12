import sys
import numpy as np
from math import *
import scipy as sp
import scipy.interpolate
from pau_parameters import *
#from gen_texp import *

#Functions.............................
def err_magnitude(m):

	#Computing photons number per pixel
	N_sig = (3631.0 * 1.51 * pow(10, 7) * np.power(10,-0.4 * m) * in_r * texp * tel_surface)/ n_pix #object photons per pixel
	N_sky = in_sky * texp * tel_surface * pix_size #sky photons per pixel

	#Computing Signal to Noise (StoN)
	Signal = N_sig * n_pix * n_exp 
	Noise =  np.sqrt(n_pix * n_exp * (N_sig + N_sky + RN * RN))	
	StoN = Signal / Noise 

	#Translating from S/N to magnitude uncertainty
	err_m_obs = 2.5 * np.log10(1 + 1 / StoN)

	#Inserting extra noise to simulate calibration errors
	noise_ctn = 2.5 * log10(1 + 0.02)
	err_m_obs = np.sqrt((err_m_obs * err_m_obs) + (noise_ctn * noise_ctn))

	rn = np.random.normal(size = n_filt)
	m_obs = m + rn * err_m_obs	

	return m_obs, err_m_obs 

def noiseless_test():

	print "\nNoiseless test:"

	f_noiseless_test = open("noiseless_test.cat", "r") 
	f_noiseless = open(f_noiseless_file, "r") 
	
	for test in f_noiseless_test:
		test = np.array(test.split(), dtype = 'string')	
		cat = np.array(f_noiseless.readline().split(), dtype = 'string')
		for i in range(1, len(test) - 2):
			if(test[i] != cat[i]):
				print "Failed"
				print "Test value %s is different from %s" % (test[i], cat[i])
				return
		
	f_noiseless_test.close()
	f_noiseless.close()

	print "Ok"
	return

def noisy_test():

	print "\nNoisy test:"

	f_noisy_test = open("noisy_test.cat", "r") 
	f_noisy = open(f_noisy_file, "r") 

	for test in f_noisy_test:
		test = np.array(test.split(), dtype = 'string')	
		cat = np.array(f_noisy.readline().split(), dtype = 'string')
		for i in range(2, len(test) - 2, 2):
			if(test[i] != cat[i]):
				print "Failed"
				print "Test value %s is different from %s" % (test[i], cat[i])
				return

	f_noisy_test.close()
	f_noisy.close()

	print "Ok"
	return

#Loading input catalog file.............
cat = np.loadtxt(cat_folder + incat_file + incat_file_extention, usecols = (index_id, index_m, index_z, index_t), unpack = 'True')

id = cat[0]
m_ref = cat[1]
z = cat[2]
t = cat[3]
n_gal = len(m_ref)

#Loading filters....................... 
filt_folder, filt_list = np.loadtxt(filt_path + filt_names_file, dtype = 'string', unpack = True)
n_filt = len(filt_list)
filt_folder = np.hstack(([ref_filt_folder], filt_folder))
filt_list = np.hstack(([ref_filt], filt_list))

#Computing in_r.........................
y_r = {}
x = {}
in_r = np.zeros(n_filt + 1)
for i in range(n_filt + 1): 
	R = np.loadtxt(filt_path + filt_folder[i] + filt_list[i] + ".res", unpack = True)

	filt = filt_list[i]

	#Defining lambda range
	x[filt] = np.arange(R[0].min(), R[0].max(), dx)

	#Interpolation
	y_r[filt] = sp.interpolate.interp1d(R[0], R[1])(x[filt]) 
	in_r[i] = (y_r[filt] / x[filt]).sum() * dx

in_r0 = in_r[0]
in_r = in_r[1:]

#Noiseless........................................................
#.................................................................

#Loading seds...........................
sed_list = np.loadtxt(sed_folder + sed_list, dtype = 'string', unpack = True)
n_sed = len(sed_list) 

#dt = (t.max() - t.min()) / (n_sed - 1)
#t_range = np.arange(t.min(), t.max() + dt, dt )
dt = 1
t_range = np.arange(0, 66)

s = {}
for sed in sed_list: s[sed] = np.loadtxt(sed_folder + sed, unpack = True)

print "# galaxies =", n_gal

f_noiseless = open(f_noiseless_file, "w") 
#Iterating over all the galaxies
for n in range(n_gal):
	sys.stdout.write("\r\tComplet: %2.2f%%" % (100 * float(n) / float(n_gal)))
	sys.stdout.flush()
	f_noiseless.write("%d " % id[n])

	s_ind_high = np.searchsorted(t_range, t[n])
	s_high = sed_list[s_ind_high]
	s_low = sed_list[s_ind_high - 1]
	coef_low = (t_range[s_ind_high] - t[n]) / dt 
	coef_high = 1 - coef_low  

	#Computing in_s................
	in_s = np.zeros(n_filt + 1)
	for i in range(n_filt + 1):
		filt = filt_list[i]
		y_s_low = sp.interpolate.interp1d(s[s_low][0] * (1 + z[n]), s[s_low][1])(x[filt]) 
		y_s_high = sp.interpolate.interp1d(s[s_high][0] * (1 + z[n]), s[s_high][1])(x[filt]) 
		y_s = coef_high * y_s_high + coef_low * y_s_low
		in_s[i] = (y_s * y_r[filt] * x[filt]).sum() * dx
		
	in_s0 = in_s[0]
	in_s = in_s[1:]

	#Generating magnitudes.........
	for i in range(n_filt):

		#Computing true magnitudes
		m = m_ref[n] + 2.5 * (log10(in_s0) + log10(in_r[i]) - log10(in_r0) - log10(in_s[i])) 
		f_noiseless.write("%8.8f " % m)

	f_noiseless.write("%4.4f %2.2f\n" % (z[n], t[n]))

f_noiseless.close()

#Noisy.........................................................
#Computing observed magnitudes: adding gaussian noise of sigma = err_m_obs
#.................................................................

#Loading exposure times................
texp = np.loadtxt(texp_file, unpack = True)

#Computing in_sky
SKY = np.loadtxt(sky_file, unpack = True) 
in_sky = np.zeros(n_filt + 1)
for i in range(n_filt + 1):
	filt = filt_list[i]
	y_sky = sp.interpolate.interp1d(SKY[0], SKY[1])(x[filt])
	in_sky[i] = (y_sky * y_r[filt]).sum() * dx
in_sky = in_sky[1:]

f_noiseless = open(f_noiseless_file, "r") 
f_noisy = open(f_noisy_file, "w") 
#Iterating over all the galaxies
for n in range(n_gal):

	cat = np.array(f_noiseless.readline().split(), dtype = 'float')
	m = np.delete(cat, (0, n_filt + 1, n_filt + 2))	
	m_obs, err_m_obs = err_magnitude(m)
	
	f_noisy.write("%d " % id[n])
	for i in range(n_filt): f_noisy.write("%4.4f %4.4f " % (m_obs[i], err_m_obs[i]))
	f_noisy.write("%4.4f %2.2f\n" % (z[n], t[n]))

f_noiseless.close()
f_noisy.close()

#Tests.................................
#noiseless_test()
#noisy_test()
