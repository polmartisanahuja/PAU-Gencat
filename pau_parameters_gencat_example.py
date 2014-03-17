from numpy import pi

#Parameters.............................
filt_path = "../../Data/Filter/"
filt_set = '_default'
#filt_set = '_redshift'
#filt_set = '_blueshift'
#filt_set = '_log'
#filt_set = '_x0_5width'
#filt_set = '_x1_5width'

filt_names_file = "PAU_filt_list" + filt_set + ".txt"

cat_folder = '../../Data/Catalog/'
#incat_file = "mock.r260.n1e6.s10.121027"
incat_file = "mock.r240.PAU.a030.s11-16.121030"
#incat_file_extention = ".txt.gz"
incat_file_extention = ".txt"
f_noiseless_file = cat_folder + incat_file + filt_set + "_noiseless.cat"
f_noisy_file = cat_folder + incat_file + filt_set + "_noisy.cat"
#f_noisy_file = cat_folder + incat_file + filt_set + "_noisy_x2texp.cat"

sed_folder = "../../Data/Template/CE_NEW/"
sed_list = "CE_NEW_MOD.txt"

#texp_file = filt_path + "PAU_x2_exp_times.txt"
texp_file = filt_path + "PAU_exp_times.txt"
sky_file = "../../Data/Template/sky_spectrum.txt"
ref_filt_folder = 'mock.r260.n1e6.s10.121027_filters/'
#ref_filt = "r_SDSS" #Reference filter
ref_filt = "i_SDSS" #Reference filter

index_id = 0 #ID of the galaxy in de catalog
index_m = 5 #This is the reference magnitude 
index_z = 1 #The true redshift
index_t = 2 #Spectral type. The associated template

#dx = 10 #Wavelength resolution in the integrals
dx = 1 #Wavelength resolution in the integrals

#Detector parameters.....................
scale = 0.265 #arcsec/pixel
RN = 5.0 #Read-out noise of the CCDS (electrons)
D_tel = 4.2 #Diametre of the telescop (meters)
aperture = 2.0 #Assumed standard size of the galaxies (arcsec2)
n_exp = 2.0 #Number of expositions
texp_tray = [45,45,50,60,75] #Exposure times (sec) in each tray of PAUCam
texp_BB = {'up':45,'g':45,'r':50,'i':75,'z':75,'y':75} #Exposure times (sec) of each of the BB
n_filt_tray = 8

tel_surface = pi * (D_tel / 2) * (D_tel / 2) #Area of the telescope's mirror 
pix_size = scale * scale #Size of the pixels in arcsec^2 
n_pix = aperture / pix_size #Numer of pixels inside the aperture 
