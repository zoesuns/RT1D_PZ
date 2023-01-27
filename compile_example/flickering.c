#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "qss.h"
#include "base.h"
#include "cosmology.h"
#include "atomic_data.h"
#include "spec.h"
#include "calc_rates.h"
#include "propagate_flickering.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_NCELLS 5000
struct stat st = {0};

int main(int argc, char *argv[]){
    init_load_recomb_data();
    init_load_cosmology_data();
    init_eVarray();

    //char outdir[]="/home/hqchen/proximity_zone/evolve_los/B40E/a=0.1008/lightray_h1803774_n0/";
    char fname[900];
    sprintf(fname,"%slos.info",argv[1]);

    /*
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    } */

    char *eptr;
    double zuni=strtod(argv[2],&eptr);
    printf("%f\n", zuni);
    int no=6;
    double tOutList[]={ 1e5*yr2s, 3e5*yr2s, 1e6*yr2s, 3e6*yr2s, 1e7*yr2s, 3e7*yr2s};

    outArrayBlock *specBlockIn;
    int nepoch=10;
    specBlockIn=output_alloc(NUM_ENERGY_BINS, nepoch); 
    specBlockIn->current=nepoch;
    for (int i=0; i<specBlockIn->current; i++){
        specBlockIn->t[i]=i*1e6*yr2s;
        new_powerlaw_spec(specBlockIn->y[i], 1e57*int(i%2), 1.5);
    }




    int tmp, nc;
    double r[MAX_NCELLS], dr[MAX_NCELLS], nH[MAX_NCELLS], nHe[MAX_NCELLS], xHI[MAX_NCELLS], xHeI[MAX_NCELLS], xHeII[MAX_NCELLS], T[MAX_NCELLS], vlos[MAX_NCELLS];
    nc=count_nrow_data_file(fname, '#', &tmp);
    //nc=100;
    read_array_from_txt(fname, nc, '#', 9, r, dr, nH, nHe, xHI, xHeI, xHeII, T, vlos);
    cellProp cell[nc];
    double bkg[nc][6];

    for (int k=0; k<nc; k++){
        cell[k].dist_pMpc=r[k];
        cell[k].dr_pMpc=dr[k];
        cell[k].nH=nH[k];
        cell[k].nHe=nHe[k];
        cell[k].xHI=xHI[k];
        cell[k].xHeI=xHeI[k];
        cell[k].xHeII=xHeII[k];
        cell[k].T=T[k];
        cell_parse_info(&cell[k]);
        printf("cell #%d\n", k);
        init_photo_bkg(&cell[k], bkg[k]);
        //printf("bkg %e %e %e %e\n", bkg[k][0], bkg[k][1], bkg[k][2], bkg[k][3]);
        //double dist=0.01*k+0.1;
        //H_only_cell(&cell[k], dist, 0.01); 
    }

    propagate_lc_los(nc, cell, specBlockIn, zuni, no, tOutList, bkg, argv[1], "test", 1, 100); // 0 means do not trim similar spectra, 100 means output cell and incidental spec every 100 cells


    free_loaded_recomb_data();
    free_loaded_cosmology_data();


    return 0;
}
