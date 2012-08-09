import filters as ft
import numpy as np
import os as os

filt_folder = "/Users/pmarti/Dropbox/Tesi/filt/filters_120_20120612/" #Folder where the computed bands will be located 
eff_filt_folder = "/Users/pmarti/Dropbox/Tesi/filt/filters_120_20120612_eff/"
trans_folder = "/Users/pmarti/Dropbox/Tesi/filt/PAU_trans_curves/"

files = os.listdir(filt_folder)

#Compute filters..........................................................
for name in files:
	lam, R = np.loadtxt(filt_folder + name, unpack = True)
	eff_R = ft.gen_effect_filt(lam, R, trans_folder)
	filt = np.array([lam, eff_R])
	np.savetxt(eff_filt_folder + name[:-4] + ".res", filt.T, fmt = ["%2.2f", "%5.5f"])

#Plot set of narrow filt.................................................
ft.plot_set_filt(eff_filt_folder)
