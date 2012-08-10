#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gennoisecat.h"

//FUNCTIONS DECLARATION..............................................................................................................
FILE *file_init(char file_name[1000]);
float mag_err(int i, float t_exp, float *in_r, float *in_cel, float mag_ab);
float rnd_gauss();
float interpol(float x1, float x0, float y1, float y0, float x);
void mock_noise(int num_filt, float *in_r, float *in_cel);
void integrals(float *in_r, float *in_cel);
void exp_times_reader(int num_filt);

float magnitude;
float *t_exp;

//Derived instrument parameters
float tel_surface; //Area del telescopi
float pix_size; //Tamany pixel en arcsec^2
float n_pix; //Número de pixels dins l'apertura
	
main()
{
	FILE *filt; //apuntador al arxiu on hi ha la resposta del filtre
	FILE *nom_filt; //apuntador al arxiu on hi ha el nom de tots els filtres

	char nom_fich_filt_a[1000];
	char nom[1000];

	int i, filt_index, num_filt = 0; 

	float *in_r, *in_cel;
	
	tel_surface = 3.145926*(D_tel/2)*(D_tel/2); //Area del telescopi
	pix_size = scale*scale; //Tamany pixel en arcsec^2
	n_pix = aperture/pix_size; //Número de pixels dins l'apertura

	//Assignant valor a num_filt......................................................................................................
	nom_filt = fopen(nom_filt_file , "r"); //apertura arxiu amb noms dels filtres
	while(fscanf(nom_filt, "%s\n", nom_fich_filt_a) != -1) num_filt++;
	fclose(nom_filt);

	//Set memory for pointers.........................................................................................................
	//Labels are filters
	in_r = (float*) calloc(num_filt, sizeof(float)); 
	in_cel = (float*) calloc(num_filt, sizeof(float));
	t_exp = (float*) calloc(num_filt, sizeof(float));
	
	printf("\n***********************************************\n");
	printf("*                                             *\n");
	printf("*     NOISY CATALOG GENERATOR by P. Marti     *\n");
	printf("*                 Version 1                   *\n");
	printf("*                                             *");
	printf("\n***********************************************\n");
	
	//Calcular la integral convolució filtres-cel......................................................................................
	printf("\n---INTEGRALS---------\n");
	integrals(in_r, in_cel);
	
	//Read exposure times..............................................................................................................
	printf("\n---EXP TIMES---------\n");
	exp_times_reader(num_filt);
	
	//Write noisecat...................................................................................................................
	printf("\n---MOCK NOISE---------\n");
	mock_noise(num_filt, in_r, in_cel);

	return 0;
}

//FUNCTIONS DEFINITION..................................................................................................................
void integrals(float *in_r, float *in_cel)
{
	/*This rutine set all the values of the in_r, in_cel, amp_filt, max_filt, lam_filt arrays*/
	
	FILE *cat, *nom_filt, *filt, *sed;
	int k, error, filt_index;
	char nom_fich_filt_a[1000], nom_fich_filt_b[1000];
	float x, *y, *y_sed; 
	float *xn, *xn1, *yn, *yn1;
	
	xn=(float*) calloc(3, sizeof(float)); /*2 parells de valors tan dels filtres (0) com del sed (1) per interpolar entre ells. 
	El 0 també es farà servir per trobar la mag limit corresponent. xn son els primers punts i xn1 els seguents*/
	xn1=(float*) calloc(3, sizeof(float));
	yn=(float*) calloc(3, sizeof(float));
	yn1=(float*) calloc(3, sizeof(float));
	y=(float*) calloc(2, sizeof(float)); //El 0 és pel filtre i el 1 pel sed.
	y_sed=(float*) calloc(2, sizeof(float)); 
	
	cat = fopen(cat_out_file , "w"); //Creem el arxiu on anirà el cataleg amb el tall en error de magnituds
	fprintf(cat,"#Fields:\n#1:ID\n");
	
	nom_filt = fopen(nom_filt_file , "r"); //apertura arxiu amb noms dels filtres

	k=0;
	filt_index = 0;
	while(fscanf(nom_filt, "%s\n", nom_fich_filt_a) != -1 ) //Bucle per integrar sobre tots els filtres
	{
		strcpy(nom_fich_filt_b,nom_filt_folder);
		strcat(nom_fich_filt_b, nom_fich_filt_a);

		fprintf(cat,"#%d:m_obs %s\n",k+2, nom_fich_filt_a);
		k++;
		fprintf(cat,"#%d:err_m_obs %s\n",k+2, nom_fich_filt_a);
		
		filt=file_init(nom_fich_filt_b); //obrim l'arxiu del filtre corresponent
		fscanf(filt, "%f %f\n", &xn[0] , &yn[0]); //lectura 2 punts
		error=fscanf(filt, "%f %f\n", &xn1[0] , &yn1[0]);
			
		sed=fopen(sky_file , "r"); 
		fscanf(sed, "%f %f\n", &xn[1] , &yn[1]);
		fscanf(sed, "%f %f\n", &xn1[1] , &yn1[1]);
		
		x=xn[0]; //x inicial
		y[0]=yn[0]; //y[0] inicial
		
		in_r[filt_index]=0; 
		in_cel[filt_index]=0;
		
		while(error != -1) //Bucle per integrar sobre tot el filtre
		{
			while(x > xn1[1]) //SED
			{
				yn[1]=yn1[1];
				xn[1]=xn1[1];
				fscanf(sed, "%f %f\n", &xn1[1] , &yn1[1]);
			} 
			y[1]=interpol(xn1[1], xn[1], yn1[1], yn[1], x);
					
			in_r[filt_index]+=y[0]*(delta/x); //integral filtre
			in_cel[filt_index]+=y[0]*y[1]*delta; //integral filtre i sed
			x=x+delta; //nou punt x
			
			while((x > xn1[0]) && (error != -1)) //Filtre
			{
				yn[0]=yn1[0];
				xn[0]=xn1[0];
				error=fscanf(filt, "%f %f\n", &xn1[0] , &yn1[0]);
			}
			if(error != -1) y[0]=interpol(xn1[0], xn[0], yn1[0], yn[0], x);
		}
		fclose(filt); 
		fclose(sed);
		
		filt_index++;
		k++;
	}
	fclose(nom_filt);
	fclose(cat);
}

void exp_times_reader(int num_filt)
{
	/*Reads exposure times from a file and then fills the t_exp[] array with them.*/	

	FILE *temps_exp;
	int i;
	
	temps_exp=fopen(temps_exp_file, "r"); 
	for(i=0; i < num_filt; i++) fscanf(temps_exp, "%f\n", &t_exp[i]);
	printf("Hello\n");	
	fclose(temps_exp);
}

void mock_noise(int num_filt, float *in_r, float *in_cel)
{
	/*This routine takes the noisless mock catalog and convert it into 
	a noiseful one. Noise is generated through a gauss distribution.
	It also set the values of the mag_lim array. Filters with a err_mag bigger than
	sigma_lim will have a magnitude value -99. and an err_mag = 0*/
	
	FILE *mock, *cat;
	int i, num_gal, N_gal = 0, n_gal = 0;
	float r, z, type, other1, other2, mag_ab, sigma_mag_ab;
	char let = 'a';
	int error, pos;
	
	printf("   Computing number of objects in mock...\n");
	
	cat = fopen(cat_out_file, "a");
	
	pos = 2 * num_filt + 1;
	
	fprintf(cat,"#%d:z_true\n#%d:m_obs prior\n#%d:template\n#%d:E(B-V)\n#%d:Mr\n", pos + 1, pos + 2, pos + 3, pos + 4, pos + 5, pos + 6 );
	
	//Count number of lines in mock
	mock = file_init(cat_in_file); 
	while(error = fscanf(mock, "%c", &let) != -1) if(let == '\n') N_gal++;
	fclose(mock);
	printf("   Total number of obj.:\t%d\n", N_gal);
	printf("   Computing observed magnitudes for all objects...\n");
	mock = file_init(cat_in_file); //Initialitzation of the noiseless cat
	while(fscanf(mock, "%d", &num_gal) != -1)
	{
		n_gal++;
		//printf("   Process completed:\t%2.2f%%\r", 100*((float)n_gal/(float)N_gal));
		//fflush(stdout);
		
		for(i=0; i < num_filt; i++)
		{
			fscanf(mock, " %f", &mag_ab);
			sigma_mag_ab = mag_err(i,t_exp[i], in_r, in_cel, mag_ab);

			mag_ab = magnitude;
			 
			if(i==0) fprintf(cat, "%d ", num_gal);
			fprintf(cat, " %4.4f %4.4f", mag_ab, sigma_mag_ab);
			if(i == num_filt-1) {
				fscanf(mock, " %f %f", &z, &type);
				fprintf(cat, " %4.4f %2.2f\n", z, type); //(FOR LEPHARE AND BPZ)			
			}
		}
	}
	fclose(cat);
	printf("\n");
}

FILE *file_init(char file_name[100])
{
	/*Returns a pointer to the file "in" which has been initialized to the first 
	row where data is. This means that all headers of the file that start with # 
	are avoided*/
	
	FILE *in;
	int N, i;
	char caracter1, caracter2;
	
	//Count number of lines begining with # (headers)
	if ((in=fopen(file_name, "r")) == NULL) printf ("Error opening file \"%s\".\n",file_name);
	N=0;
	fscanf(in, "%c", &caracter1);
	while(caracter1=='#') {
		N++;
		fscanf(in, "%c", &caracter2);
		while(caracter2!='\n') fscanf(in, "%c", &caracter2);
		fscanf(in, "%c", &caracter1);
	}
	fclose(in);
	
	//Read rows until data will be found
	if ((in=fopen(file_name, "r")) == NULL) printf ("Error opening file %s.\n",file_name);
	i=0;
	while(i != N) {
		fscanf(in, "%c", &caracter1);
		if(caracter1 == '\n') i++;
	}
		
	return in;
}

float mag_err(int i, float t_exp, float *in_r, float *in_cel, float mag_ab)
{
	//Computes the error of the magnitude for a given magnitude, 
	//the exposure time in some band, the integral of the band throughput
	//and the inegral of the band throughput multiplied by the sky brigthness

	float N_sig, N_cel, StoN;
	float sigma_mag_ab, noise_ctn;
	float r;			

	//Computing photons number
	N_sig = (3631.0*1.51*pow(10, 7)*pow(10,-0.4*mag_ab)*in_r[i]*t_exp*tel_surface)/n_pix; //object's photons per pixel
	N_cel = in_cel[i]*t_exp*tel_surface*pix_size; //sky's photons per pixel
	
	//Computing Signal to Noise (StoN)
	StoN = (sqrt(n_pix*n_exp)*N_sig)/sqrt(N_sig+N_cel+R*R);
	sigma_mag_ab = 2.5*log10(1+1/StoN);
	noise_ctn = 2.5*log10(1+noise);
	sigma_mag_ab = sqrt((sigma_mag_ab*sigma_mag_ab)+(noise_ctn*noise_ctn)); //afegim un soroll uniforme a tots els filtres

	//Add noise
	r = rnd_gauss(); //r: random value from a normal distribution
	magnitude = mag_ab + r * sigma_mag_ab;

	return sigma_mag_ab;
}

float rnd_gauss()
{
	/*Returns a random value from a normal distribution*/

	float r1, r2, r;
	
	r1 = (float)rand()/RAND_MAX;//Generador destribucio gaussiana
	r2 = (float)rand()/RAND_MAX;
	r = sqrt(-2*log(r1))*cos(2*3.145926*r2);
	r = r*preci;
	r = (ceil(r))/preci;

	return r;
}

float interpol(float x1, float x0, float y1, float y0, float x)
{
	/*Returns the linearly interpolated y between two points in the plane*/

	float m, n;
	
	m=(y1-y0)/(x1-x0);
	n=y0-m*x0;
	return m*x+n;
}
