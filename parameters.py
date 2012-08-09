#Parameters.............................
incat_file = "mock.r255.PAU.a010.s12.110521_reduced"
incat_file_extention = ".dat"
f_noiseless_file = incat_file + "_noiseless.cat"

sed_folder = "./CE_NEW/"
sed_list = "CE_NEW_MOD.txt"

filt_folder = "./PAU_filt/"
filt_list = "PAU_filt_list.txt"

ref_filt = "r_SDSS.res" #Reference filter

index_id = 0 #ID of the galaxy in de catalog
index_m = 5 #This is the reference magnitude 
index_z = 1 #The true redshift
index_t = 2 #Spectral type. The associated template

dx = 0.1 #Wavelength resolution in the integrals
