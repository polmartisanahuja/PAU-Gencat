#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define delta 0.1 //El pas d'integració
#include "gennoislesscat.h"

FILE* file_init(char file_name[1000]);
float interpol(float x1, float x0, float y1, float y0, float x);
void set_type(float type, float *coef, char *sed1, char *sed2);
void filt_path(char *filt_b, char *filt_a);
void inter_sed_point(FILE *sed, int i, float *x1, float *x0, float *y1, float *y0, float x, float z);
int inter_filt_point(FILE *filt, int *error, int i, float *x1, float *x0, float *y1, float *y0, float x);

main()
{
	FILE *mock; //Input catalog
	FILE *cat; //Output catalog
	FILE *filt; //Filter response
	FILE *nom_filt; //Filter names
	FILE *nom_sed; //Sed names
	FILE *sed1; //Spectral template. It must cover all filters range
	FILE *sed2; //Spectral template. It must cover all filters range  
	
	char nom_fich_filt_a[1000];
	char nom_fich_filt_b[1000];
	char nom_fich_sed1[1000];
	char nom_fich_sed2[1000];
	char line[1000];

	int i, N_gal = 0, *error, filt_index, sed_index1, sed_index2, num_gal = 0, num_filt = 0; //index filtr, sed & obj
	float type; //Real spectral index
	float *in_r, *in_sed; //in: integral r:filter sed:template
	float x, *y, *y_sed, *xn, *xn1, *yn, *yn1; 
	float id, ebv, Mr, z, mi, mag_ab;
	float *coef;
	float intro[5], sdss[2], des[6], pau[40], F814, mass[3];

	xn = (float*) calloc(3, sizeof(float)); //0:filter 1:sed1 2:sed2
	xn1 = (float*) calloc(3, sizeof(float)); //0:filter 1:sed1 2:sed2
	yn = (float*) calloc(3, sizeof(float)); //0:filter 1:sed1 2:sed2
	yn1 = (float*) calloc(3, sizeof(float)); //0:filter 1:sed1 2:sed2
	y = (float*) calloc(2, sizeof(float)); //0:filter 1:sed
	y_sed = (float*) calloc(2, sizeof(float)); //0:sed1 1:sed2
	coef = (float*) calloc(2, sizeof(float)); //0:sed1 1:sed2
	error = (int*) calloc(1, sizeof(int));

	printf("\n***********************************************\n");
	printf("*                                             *\n");
	printf("*   NOISELESS CATALOG GENERATOR by P. Marti   *\n");
	printf("*                 Version 1                   *\n");
	printf("*                                             *");
	printf("\n***********************************************\n");

	//Computing number of filters....................................................................
	nom_filt = fopen(nom_filt_file, "r"); //apertura arxiu amb noms dels filtres
	while(*error = fscanf(nom_filt, "%s\n", nom_fich_filt_a) != -1) if(*error !=-1) num_filt++;
	fclose(nom_filt);
	printf("\n\tTotal Number of filters = %d\n",num_filt);
		
	//Openinf output catalog.........................................................................
	cat = fopen(cat_out_file , "w"); 
	fprintf(cat,"#1:ID\n#2:z_true\n#3:t_true\n#4:E(B-V)\n#5:Mr\n");
	
	//FILTERS MODULE.................................................................................
	nom_filt = fopen(nom_filt_file, "r"); //apertura arxiu amb noms dels filtres

	//Loop over filters..............................................................................
	in_r = (float*) calloc(num_filt + 1, sizeof(float));
	i = 0;
	do{ 
		if(i == 0) filt_path(nom_fich_filt_b, nom_filt_seed);
		else filt_path(nom_fich_filt_b, nom_fich_filt_a);
		
		//filt_path(nom_fich_filt_b, nom_fich_filt_a);
		//printf("%s\n", nom_fich_filt_a);
		
		if(i != 0) fprintf(cat,"#%d:%s\n",i + 5, nom_fich_filt_a);
		
		//filt = fopen(nom_fich_filt_b , "r"); //obrim l'arxiu del filtre corresponent
		filt = file_init(nom_fich_filt_b); //obrim l'arxiu del filtre corresponent

		fscanf(filt, "%f %f\n", &xn[0] , &yn[0]); //lectura 2 punts
		*error = fscanf(filt, "%f %f\n", &xn1[0] , &yn1[0]);
		x = xn[0]; //x inicial
		y[0] = yn[0]; //y[0] inicial
		
		in_r[i] = 0;
		while(*error != -1){ //Bucle per integrar sobre tot el filtre
			in_r[i] += y[0]*(delta/x); //integral filtre
			x+=delta; //nou punt x
			
			*error = inter_filt_point(filt, error, 0, xn1, xn, yn1, yn, x);
			if(*error != -1) y[0] = interpol(xn1[0], xn[0], yn1[0], yn[0], x);
		}
		fclose(filt); //tancament arxiu filtre
		i++;
	}while(fscanf(nom_filt, "%s\n", nom_fich_filt_a) != -1);

	fclose(nom_filt);
	
	//CATALOG MODULE.................................................................................
	//Count number of objects in the file
	mock = file_init(cat_in_file);
	while(fgets(line, 1000, mock) != NULL) N_gal++;
	fclose(mock);
	printf("\tTotal number of objects = %d\n", N_gal);

	mock = file_init(cat_in_file);
	
	//Loop over objects of mock catalog..............................................................
	in_sed = (float*) calloc(num_filt + 1, sizeof(float));
	//while(fscanf(mock, "%f %f %f %f %f %f\n", &id, &mi, &z, &type, &ebv, &Mr) != -1){
	while(fscanf(mock,"%f ", &intro[0]) != -1){
		num_gal++;
		
		for(i=1; i<5; i++) fscanf(mock,"%f ",&intro[i]);
		fscanf(mock,"%f ", &sdss[0]);
		for(i=0; i<6; i++) fscanf(mock,"%f ",&des[i]);
		for(i=0; i<40; i++) fscanf(mock,"%f ",&pau[i]);
		fscanf(mock,"%f ", &F814);
		fscanf(mock,"%f ", &sdss[1]);
		for(i=0; i<2; i++) fscanf(mock,"%f ",&mass[i]);
		fscanf(mock,"%f\n",&mass[2]);
		
		//printf("id = %0.0f\n", intro[0]);
		id = intro[0];
		z = intro[1];
		type = intro[2];
		ebv = intro[3];
		Mr = intro[4];
		mi = sdss[1];
		
		//printf("%f\n", mi);
		
		//Calcular la integral convolució filtres-SED i filtres sols.................................
		nom_filt = fopen(nom_filt_file, "r"); //apertura arxiu amb noms dels filtres
			
		set_type(type, coef, nom_fich_sed1 , nom_fich_sed2);
		//printf("%f %f\n", coef[0], coef[1]); 
		//Loop over filters..........................................................................
		i = 0;
		do{ 
			if(i == 0) filt_path(nom_fich_filt_b, nom_filt_seed);
			else filt_path(nom_fich_filt_b, nom_fich_filt_a);
			
			sed1 = fopen(nom_fich_sed1 , "r");
			sed2 = fopen(nom_fich_sed2 , "r");
			//printf("%s %s\n", nom_fich_sed1, nom_fich_sed2);
			xn1[1] = 0, xn1[2] = 0; //Forcing SEDs to start at wavelenght 0 
			
			filt = file_init(nom_fich_filt_b); //obrim l'arxiu del filtre corresponent
			//printf("%s\n", nom_fich_filt_b);
			fscanf(filt, "%f %f\n", &xn[0] , &yn[0]); //lectura 2 punts
			*error = fscanf(filt, "%f %f\n", &xn1[0] , &yn1[0]);
			x = xn[0]; //x inicial
			y[0] = yn[0]; //y[0] inicial
			
			//Loop over all the filter range.........................................................
			in_sed[i] = 0;
			while(*error != -1){ 
				//SED1
				inter_sed_point(sed1, 1, xn1, xn, yn1, yn, x, z);
				y_sed[0] = interpol(xn1[1], xn[1], yn1[1], yn[1], x);
				
				//SED2
				inter_sed_point(sed2, 2, xn1, xn, yn1, yn, x, z); 
				y_sed[1] = interpol(xn1[2], xn[2], yn1[2], yn[2], x);
		
				y[1] = coef[0]*y_sed[0] + coef[1]*y_sed[1]; //Combinació lineal seds
				in_sed[i] += y[0]*y[1]*(x*delta); //integral filtre i sed
				x += delta; //nou punt x
				
				//FILTER
				*error = inter_filt_point(filt, error, 0, xn1, xn, yn1, yn, x);
				if(*error != -1) y[0] = interpol(xn1[0], xn[0], yn1[0], yn[0], x);
				
			}
			fclose(filt), fclose(sed1), fclose(sed2);
			//printf("%f ", in_sed[i]);	
			i++;
		}while(fscanf(nom_filt, "%s\n", nom_fich_filt_a) != -1);
		
		fclose(nom_filt);

		//Printing catalog with all magnitudes in different filters....................................
		fprintf(cat, "%d ", (int)id);
		printf("\tProcess completed:\t%2.2f%%\r", 100*((float)num_gal/(float)N_gal));
		fflush(stdout);

		for(i=1; i <= num_filt; i++){
			mag_ab = mi+2.5*(log10(in_sed[0])+log10(in_r[i])-log10(in_r[0])-log10(in_sed[i]));
			fprintf(cat, "%4.4f ", mag_ab);
		}
		fprintf(cat, "%4.4f %2.2f\n", z, type);
	}
	fclose(cat);
	fclose(mock);
	
	printf("\n");

	return 0;
}

FILE *file_init(char file_name[1000])
{
	/*Returns a pointer to the file "in" which has been initialized to the first 
	row where data is. This means that all headers of the file that start with # 
	are avoided*/
	
	FILE *in;
	int N = 0, i = 0;
	char caracter1, caracter2;
	
	//Count number of lines begining with # (headers)
	if ((in=fopen(file_name, "r")) == NULL) printf ("Error opening file \"%s\".\n",file_name);

	fscanf(in, "%c", &caracter1);
	while(caracter1=='#') {
		N++;
		fscanf(in, "%c", &caracter2);
		while(caracter2!='\n') fscanf(in, "%c", &caracter2);
		fscanf(in, "%c", &caracter1);}
	fclose(in);
	
	//Read rows until data will be found
	if ((in=fopen(file_name, "r")) == NULL) printf ("Error opening file %s.\n",file_name);

	while(i != N) {
		fscanf(in, "%c", &caracter1);
		if(caracter1 == '\n') i++;}
		
	return in;
}

void filt_path(char *filt_b, char *filt_a)
{
	/*Append filt_a behind filt_b*/
	
	strcpy(filt_b,nom_filt_folder);
	strcat(filt_b, filt_a);
	
	//printf("")
}
	
void set_type(float type, float *coef, char *sed1, char *sed2)
{
	/*Gives values to the "coef" variable acording to "type" variable. It also set paths for
	both templates that represent the correspondent spectral type.*/

	FILE *nom_sed;
	char nom_fich_sed1_a[1000];
	char nom_fich_sed2_a[1000];
	int sed_index1, sed_index2;
	
	nom_sed = fopen(nom_sed_file , "r"); //apertura arxiu amb noms dels seds
	sed_index1 = -1;
	while(sed_index1 != (int)type){
		fscanf(nom_sed, "%s\n", nom_fich_sed1_a); //lectura index i nom del primer sed
		strcpy(sed1,nom_sed_folder);
		strcat(sed1, nom_fich_sed1_a);
		sed_index1+=1;
		coef[0]=(float)1-type+sed_index1;
	}
	if(sed_index1 != 65) fscanf(nom_sed, "%s\n", nom_fich_sed2_a); //lectura index i nom del primer sed. Cal arreglar el 65 que es el nombre de seds del FJC.
	else strcpy(nom_fich_sed2_a, nom_fich_sed1_a);
	strcpy(sed2,nom_sed_folder);
	strcat(sed2, nom_fich_sed2_a);
	coef[1]=(float)type-sed_index1;
	fclose(nom_sed);
	
	//printf("sed_index1=%d\n", sed_index1);
	//printf("sed1=%s\n", nom_fich_sed1_a);
	//printf("sed2=%s\n", nom_fich_sed2_a);
	//printf("coef[0]=%f coef[1]=%f\n", coef[0], coef[1]);
}

void inter_sed_point(FILE *sed, int i, float *x1, float *x0, float *y1, float *y0, float x, float z)
{
	/*Finds two points in the "sed" file that enclose the point x, then changes pointers values x & y for
	both 0 & 1.*/
	
	while(x > x1[i]){ 
		y0[i] = y1[i], x0[i] = x1[i];
		fscanf(sed, "%f %f\n", &x1[i] , &y1[i]);
		x1[i] = x1[i]*(1+z); //Wavelength Redshifted
	}
}

int inter_filt_point(FILE *filt, int *error, int i, float *x1, float *x0, float *y1, float *y0, float x)
{
	/*Finds two points in the "filt" file that enclose the point x, then changes pointers values x & y for
	both 0 & 1. If the end of filter is found, then "error" pointer tooks value -1 and loop is
	aborted*/
	
	while((x > x1[i]) && (*error != -1)){ 
		y0[i] = y1[i], x0[i] = x1[i];
		*error = fscanf(filt, "%f %f\n", &x1[i] , &y1[i]);
	}
	return *error;
}

float interpol(float x1, float x0, float y1, float y0, float x)
{
	/*Returns the linearly interpolated value y between two points in the plane (x, y)*/

	float m, n;
	
	m=(y1-y0)/(x1-x0);
	n=y0-m*x0;
	return m*x+n;
}

