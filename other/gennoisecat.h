//Physical parameter
//Color noisy cat generator
#define noise 0.02 // N/S de soroll constant que s'afegeix al error de les magnituds
#define scale 0.27 // arcsec/pixel del detector
#define R 5.0 //soroll de lectura dels ccds (electrons)
#define D_tel 4.2 //Diametre telescopi (meters)
#define aperture 2.0 //Tamany galaxies. De moment considerarem que totes tenen el mateix tamany (arcsec2)
#define n_exp 2.0 //Número de exposicions
#define filter_sys 42

//Numerical parameters
#define delta 0.1 //El pas d'integracio
#define preci 100000 //Dona el número de decimals que tindran els números aleatoris del 0 al 1 

//File paths
char cat_in_file[1000] = "mock.r255.PAU.a010.s12.110521_reduced_noiseless_c.dat";
char cat_out_file[1000] = "mock.r255.PAU.a010.s12.110521_reduced_noisy_c.dat";
char nom_filt_file[1000] = "./PAU_filt/PAU_filt_list.txt";
char nom_filt_folder[1000] = "./PAU_filt/";
char sky_file[1000] = "sky_spectrum.txt";
char temps_exp_file[1000] = "exp_times.txt";
