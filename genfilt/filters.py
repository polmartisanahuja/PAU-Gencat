import numpy as np
import scipy as sp
import scipy.interpolate
import os as os
import sys
import matplotlib.pyplot as plt

def gen_filt(x0, Dx, dx, x_step):
	"""Generates a single rectangle filter with main body 
	width Dx and wings width dx.The main body starts at x0 """
	
	if(Dx == 0): sys.exit("Filter width must be non-zero")
	if(dx == 0): dx = 0.01 * Dx 
	
	y = np.array([0,0.95,0.95,0])
	x = np.array([x0 - dx, x0, x0 + Dx, x0 + Dx + dx])
	
	new_x = np.arange(x[0],x[3] + x_step, x_step)
	new_y = interpol(x, new_x, y)
	
	return new_x, new_y
	
def gen_effect_filt(lam, R, trans_curves_folder):
	"""Returns the filter x,y convoluted with all the transmission curves in
	the ~/routines/trans_curves folder"""
	
	files = os.listdir(trans_curves_folder)
	
	eff_R = R
	for name in files:
		tx, ty = np.loadtxt(trans_curves_folder + name, unpack = True)
		t = sp.interpolate.interp1d(tx,ty)
		eff_R = eff_R * t(lam)
	
	return eff_R
	
def interpol(x, newx, y):
	"""Returnes the array y once has been interpolated into the newx array. 
	The interpolation is made linearly"""
	
	f = sp.interpolate.interp1d(x,y)
	
	return f(newx)
	
def plot_set_filt(folder_name):
	"""Plot the filters are in the folder_name"""

	try: os.remove(folder_name + ".DS_Store")
	except: pass
	
	files = os.listdir(folder_name)

	for name in files:
		lam, R = np.loadtxt(folder_name + name, unpack = True)
		plt.fill(lam, R, alpha = 0.75)
		plt.xlabel("$\lambda$ (A)")
		plt.ylabel("Throughput")
		plt.savefig(folder_name + "plot_filt.png")
		
def plot_set_trans_curves(folder_name):

	files = os.listdir(folder_name)

	try: os.remove(folder_name + ".DS_Store")
	except: pass
	
	
	for name in files:
		lam, T = np.loadtxt(folder_name + name, unpack = True)
		plt.plot(lam, T)
		plt.xlabel("$\lambda$ (A)")
		plt.ylabel("Transmission")

		
def gen_set_effect_filt(folder_name, trans_curves_folder):
	
	try: os.remove(folder_name + ".DS_Store")
	except: pass
	
	files = os.listdir(folder_name)
	folder_name_effective = folder_name[:-1] + "_effective"
	os.system("mkdir " + folder_name_effective)
	
	for name in files:
		lam, R = np.loadtxt(folder_name + "/" + name, unpack = True)
		eff_R = gen_effect_filt(lam, R, trans_curves_folder)
		filt = np.array([lam, eff_R])
		np.savetxt(folder_name_effective + "/" + name, filt.T, fmt = ["%2.2f", "%5.5f"])
	
	
