import filters as ft
import numpy as np
import os as os

filt_folder = "/Users/pmarti/Professional/Studies/3_PhD/Tesi/filt/PAU_wide/" #Folder where the computed bands will be located 
trans_folder = "/Users/pmarti/Professional/Studies/3_PhD/Tesi/filt/PAU_trans_curves/"

#Create folder...........................................................
try:
	os.mkdir(filt_folder)
except:
	print "Folder already exists:" + filt_folder
	
#Compute effective filters...............................................
ft.gen_set_effect_filt(filt_folder, trans_folder)

#Plot set of narrow filt.................................................
ft.plot_set_filt(filt_folder[:-1] + "_effective/")
	