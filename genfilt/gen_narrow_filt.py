import filters as ft
import numpy as np
import os as os

n_filt = 42 #Number of narrow bands

x0 = 4400 #Lambda min (without wings) of the bluest filter (Amstrongs)
Dx = 100 #Width of the filters without wings (Amstrongs)
dx = 25 #Width of the filter wings (Amstrongs)
x_step = 1 #Step of the curve (Amstrongs)

filt_folder = "/Users/pmarti/Professional/Studies/3_PhD/Tesi/filt/" + "PAU_narrow_x0" + str(x0) + "_Dx" + str(Dx) + "_dx" + str(dx) + "/" #Folder where the computed bands will be located 
trans_folder = "/Users/pmarti/Professional/Studies/3_PhD/Tesi/filt/PAU_trans_curves/"

#Create folder............................................................
try:
	os.mkdir(filt_folder)
except:
	print "Folder already exists:" + filt_folder
	
#Compute filters..........................................................
for i in range(n_filt):
	x_init = x0 + i * Dx
	lam, R = ft.gen_filt(x_init, Dx , dx, x_step)
	eff_R = ft.gen_effect_filt(lam, R, trans_folder)
	filt = np.array([lam, eff_R])
	np.savetxt(filt_folder + str(x_init) + ".res", filt.T, fmt = ["%2.2f", "%5.5f"])

#Plot set of narrow filt.................................................
ft.plot_set_filt(filt_folder)
	