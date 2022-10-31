#ifndef _INCL_BASE_
#define _INCL_BASE_

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define NUM_ENERGY_BINS 80
#define NU_MIN 13.6
#define NU_MAX 1000.0

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define Mpc2cm (3.086e24)
#define eV2erg (1.60218e-12)
#define k_Boltzmann (1.380649e-16)
#define protonMass (1.6726219e-24)
#define yr2s (3.154e+7)

#define VALUE_ERROR -1

double eVarrayLeft[NUM_ENERGY_BINS];
double eVarrayRight[NUM_ENERGY_BINS];
double eVarrayCen[NUM_ENERGY_BINS];
double eVarrayWidth[NUM_ENERGY_BINS];

int init_eVarray();



int print_spec_related(double *eV, int ifWrite, char *fileName);


typedef struct cell_t{
    //properties of a cell
    double dist_pMpc,dr_pMpc;
    double dist_cm, dr_cm;

    double nH, nHe;
    double xHI, xHeI, xHeII;
    double T;

    double xHII, xHeIII;
    double nHI, nHII;
    double nHeI, nHeII,nHeIII;
    double ne, ntot;
}cellProp;

int copy_cell(cellProp *copy, cellProp original);

int cell_parse_info(cellProp *cell);

void print_cell_properties(cellProp *cell);

int write2file_ty(double *t, double **y, int num_var, int length, char *fileName);
int write_result(double d[], double y[][4], int num_var, int length, char *fileName);
 
#define MAXROW      2000
#define MAXCOL      10
 
int count_nrow_data_file(char *filename, char headerMark, int *nh);        
            
int read_1D_data_file(char *filename, int nrow, double *y, char headerMark, int nh);
int read_2D_data_file(char *filename, int nrow, double *x, double *y, char headerMark, int nh);


int write_1D_array_to_txt(double a1[], int nrow, char *filename, char* format, char* header);
int read_1D_array_from_txt(char *filename, int nrow, double a1[], char headerMark);
int write_2D_array_to_txt(double a1[], double a2[], int nrow, char *filename, char* fmt, char* header);
int read_2D_array_from_txt(char *filename, int nrow, double a1[], double a2[], char headerMark);

int write_1_4(double d[], double y[][4], int length, char *fileName);

int read_array_from_txt(char *filename, int nrow, char headerMark, int n_args, ...);

#endif
