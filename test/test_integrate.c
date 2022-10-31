#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "../base.h"
#include "../atomic_data.h"
#include "../spec.h"
#include "../calc_rates.h"
#include "../qss.h"
#include "../integrate.h"
#include "fake_cells.h"


#include <fenv.h>
#include <signal.h>

struct stat st = {0};

double err[4]={1e-2,1e-2,1e-2,1e-2};

int testing(cellProp *cell, char ouputFileName[], double *spec, char outSpecFileName[]){

    qss_system *qs;
    qs=qss_alloc(4, rate_func, NULL);


    double y[4]={cell->xHI, cell->xHeI, cell->xHeII, cell->T};
    printf("y_0=%10.4e,%10.4e,%10.4e,%10.4e\n",y[0],y[1],y[2],y[3]);

    double bkg[6]={0,0,0,0,0,0};
    
    rateParams  inParams={spec, *cell, bkg, 6};

    int iout_prev=0;

    // tmpOut does not include the initial y0
    outArrayBlock *tmpOut=output_alloc(4,100);

    // specBlockOut include the transmitted spectrum at t=1, thus specBlockOut->current should be tmpOut->current+1
    outArrayBlock *specBlockOut=output_alloc(NUM_ENERGY_BINS,100);
    specBlockOut->current=1;
    specBlockOut->t[0]=0;
    calc_transmitted_spec(specBlockOut->y[0], spec, y, cell);

    int numOut=8;
    double tOutList[]={1e5*yr2s, 3e5*yr2s, 1e6*yr2s, 3e6*yr2s, 1e7*yr2s, 3e7*yr2s, 6e7*yr2s, 1e8*yr2s};
    //int numOut=1;
    //double tOutList[]={ 1e8*yr2s};

    iout_prev=tmpOut->current;
    qss_solve_save(qs,0,tOutList[0],y,err, &inParams, tmpOut);
    save_all_spec(specBlockOut, iout_prev, spec, tmpOut, &inParams.cell_pass);

    for (int j=0; j<numOut-1; j++){
        printf("%d\n", j);
        update_inParams(&inParams, spec, y, bkg, 6);
        iout_prev=tmpOut->current;
        qss_solve_save(qs,tOutList[j],tOutList[j+1],y,err, &inParams, tmpOut);
        printf(" tmpOut->current+1= %lu\n",  tmpOut->current+1);
        save_all_spec(specBlockOut, iout_prev, spec, tmpOut, &inParams.cell_pass);
        printf("y_1=%10.4e,%10.4e,%10.4e,%10.4e\n",y[0],y[1],y[2],y[3]);
        printf("spec[0]=%e, spec[1]=%e\n", spec[0], spec[1]);
    }
    
    write2file_ty(tmpOut->t, tmpOut->y, 4, tmpOut->current-1, ouputFileName);
    write2file_ty(specBlockOut->t, specBlockOut->y, NUM_ENERGY_BINS, specBlockOut->current, outSpecFileName);
    qss_free(qs);

    return 0;
}



int main(){
    //feenableexcept(FE_INVALID | FE_OVERFLOW);

    char outdir[]="integrate_test/";

    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }


    init_eVarray();    
    init_load_recomb_data();

    double *spec;
    spec=init_spec();


    cellProp cell;

    /*
    void (*cellFunc[4]) (cellProp *cell, double dist_pMpc, double dr_pMpc);
    cellFunc[0]=cold_dense_neutral_cell;
    cellFunc[1]=hot_diffuse_ionized_cell;
    cellFunc[2]=H_only_cell;
    cellFunc[3]=He_only_cell;

    char outputNames[][100]={"cold.out", "hot.out", "H_only.out", "He_only.out"};
    char outputSpecNames[][100]={"cold.outspec", "hot.outspec", "H_only.outspec", "He_only.outspec"};
    char outn[200], outs[200];
    
    char inputSpecName[]={"spectra_test/Ndot_1e57_alpha_3.0_NHI_1e+19_NHeI_1e+16_NHeII_1e+15.dat"};
    for (int i=2; i<3; i++){
        cellFunc[i](&cell, 0.1, 0.01);
        strcpy(outn, outdir);
        strcpy(outs, outdir);
        strcat(outn, outputNames[i]);
        strcat(outs, outputSpecNames[i]);
        testing(&cell, outn, inputSpecName, outs);
    }
    */



    cold_dense_neutral_cell(&cell, 1, 0.01);
    double nHList[7]={1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1};
    char inputSpecName[][200]={"spectra_test/Ndot_1e57_alpha_1.5_NHI_1e+15_NHeI_1e+15_NHeII_1e+15.dat",
                               "spectra_test/Ndot_1e57_alpha_1.5_NHI_1e+16_NHeI_1e+16_NHeII_1e+17.dat",
                               "spectra_test/Ndot_1e57_alpha_1.5_NHI_1e+17_NHeI_1e+17_NHeII_1e+19.dat",
                               "spectra_test/Ndot_1e57_alpha_1.5_NHI_1e+18_NHeI_1e+18_NHeII_1e+19.dat",
                               "spectra_test/Ndot_1e57_alpha_1.5_NHI_1e+19_NHeI_1e+18_NHeII_1e+19.dat",
                               "spectra_test/Ndot_1e57_alpha_1.5_NHI_1e+20_NHeI_1e+18_NHeII_1e+19.dat"};


    char outname[200];
    char outspecname[200];


    for (int i=0; i<7; i++){
        for (int j=0; j<6; j++){
            read_1D_array_from_txt(inputSpecName[j], NUM_ENERGY_BINS, spec, '#');
            for (int k=0; k<NUM_ENERGY_BINS; k++){
                printf("%e ", spec[k]);
            }
            printf("\n");
            cell.nH=nHList[i];
            cell.nHe=0.078*nHList[i];
            cell_parse_info(&cell);
            sprintf(outname,"integrate_test/nH%3.1e_spec%d_8outputs.out",cell.nH,j);
            sprintf(outspecname,"integrate_test/nH%3.1e_spec%d_8outputs.outspec",cell.nH,j);
            
            testing(&cell, outname, spec, outspecname);
        }

    }

    free_loaded_recomb_data();
    return 0;
}
