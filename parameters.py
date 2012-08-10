from numpy import pi

#Parameters.............................
incat_file = "mock.r255.PAU.a010.s12.110521_reduced"
incat_file_extention = ".dat"
f_noiseless_file = incat_file + "_noiseless.cat"
f_noisely_file = incat_file + "_noisely.cat"

sed_folder = "./CE_NEW/"
sed_list = "CE_NEW_MOD.txt"

filt_folder = "./PAU_filt/"
filt_list = "PAU_filt_list.txt"

exp_t_file = "exp_times.txt"
sky_file = "sky_spectrum.txt"
ref_filt = "i_SDSS.res" #Reference filter

index_id = 0 #ID of the galaxy in de catalog
index_m = 53 #This is the reference magnitude 
index_z = 1 #The true redshift
index_t = 2 #Spectral type. The associated template

dx = 0.1 #Wavelength resolution in the integrals

#Detector parameters.....................
scale = 0.27 #arcsec/pixel
RN = 5.0 #Read-out noise of the CCDS (electrons)
D_tel = 4.2 #Diametre of the telescop (meters)
aperture = 2.0 #Assumed standard size of the galaxies (arcsec2)
n_exp = 2.0 #Number of expositions

tel_surface = pi * (D_tel / 2) * (D_tel / 2) #Area of the telescope's mirror 
pix_size = scale * scale #Size of the pixels in arcsec^2 
n_pix = aperture / pix_size #Numer of pixels inside the aperture 
