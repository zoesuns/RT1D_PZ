#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "../base.h"
#include "../cosmology.h"
#include "../atomic_data.h"
#include "../spec.h"
#include "../calc_rates.h"
#include "../propagate.h"
#include "fake_cells.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_NCELLS 3000
struct stat st = {0};

int main(int argc, char *argv[]){
    init_load_recomb_data();
    init_load_cosmology_data();
    init_eVarray();

    //char outdir[]="/home/hqchen/proximity_zone/evolve_los/B40E/a=0.1008/lightray_h1803774_n0/";
    char fname[200];
    sprintf(fname,"%slos.info",argv[1]);
    
    /*
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    } */

    char *eptr;
    double zuni=strtod(argv[2],&eptr);
    printf("%f\n", zuni);
    int no=8;
    double tOutList[]={ 1e5*yr2s, 3e5*yr2s, 1e6*yr2s, 3e6*yr2s, 1e7*yr2s, 3e7*yr2s, 6e7*yr2s, 1e8*yr2s};
    double *spec;
    spec=init_spec();
    new_powerlaw_spec(spec, 1e57, 1.5);


    int tmp, nc;
    double r[MAX_NCELLS], dr[MAX_NCELLS], nH[MAX_NCELLS], nHe[MAX_NCELLS], xHI[MAX_NCELLS], xHeI[MAX_NCELLS], xHeII[MAX_NCELLS], T[MAX_NCELLS];
    nc=count_nrow_data_file(fname, '#', &tmp);
    read_array_from_txt(fname, nc, '#', 8, r, dr, nH, nHe, xHI, xHeI, xHeII, T);
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
        printf("bkg %e %e %e %e\n", bkg[k][0], bkg[k][1], bkg[k][2], bkg[k][3]);
        //double dist=0.01*k+0.1;
        //H_only_cell(&cell[k], dist, 0.01); 
    }
    propagate_los(nc, cell, spec, zuni, no, tOutList, bkg, argv[1], "test");





    free_loaded_recomb_data();
    free_loaded_cosmology_data();


    /*
    init_load_recomb_data();
    double tmp=recomb_rate(1.0e+01 , HI);
    printf("rate=%10.4e\n",tmp);

    cellProp cell={1, 0.01, 0,0,
                    1e-4,1e-5,
                    0.6,0.7,0.2,1e3,
                    0,0, 0,0, 0,0,0,0};
    cell_parse_info(&cell);
    cellProp anotherCell=cell;
    printf("cell.nH=%10.4e\n", cell.nH);
    printf("anotherCell.nH=%10.4e\n", anotherCell.nH);
    anotherCell.nH=1000;
    printf("after modify anotherCell\n");
    printf("cell.nH=%10.4e\n", cell.nH);
    printf("anotherCell.nH=%10.4e\n", anotherCell.nH);
    

    double rec=recomb_cooling(cell.T, cell.ne, cell.nHII, cell.nHeII, cell.nHeIII);
    printf("recomb_cooling=%e\n", rec);
    double colion=colli_ioniz_cooling(cell.T, cell.ne, cell.nHI, cell.nHeI, cell.nHeII);
    printf("colli_ion_cooling=%e\n", colion);
    double colex=colli_excit_cooling(cell.T, cell.ne, cell.nHI, cell.nHeI, cell.nHeII);
    printf("colli_ex_cooling=%e\n", colex);
    double brem=brem_cooling(cell.T, cell.ne, cell.nHII, cell.nHeII, cell.nHeIII);
    printf("brem_cooling=%e\n", brem);
    double inv=inv_comp_cooling(cell.T, cell.ne, 6);
    printf("invC_cooling=%e\n", inv);

 
    double recomb=recomb_cooling(cell.T, cell.ne, cell.nHII, cell.nHeII, cell.nHeIII);
    printf("recomb=%e\n",recomb);
    calc_Cooling(&cell,6);

    free_loaded_recomb_data();
    

    init_load_recomb_data();
    double tmp=recomb_rate(1.0e+01 , HI);
    printf("rate=%10.4e\n",tmp);
    free_loaded_recomb_data();

    init_eVarray();
    print_spec_related(eVarrayWidth, 0, "ff");

    double *spec;
    spec=init_spec();

    new_powerlaw_spec(spec, 1e57, 1.5);

    print_spec_related(spec,0,"0");

    cellProp cell={1, 0.01, 0,0,
                    1e-4,1e-5,
                    0.6,0.7,0.2,1e3,
                    0,0, 0,0, 0,0,0,0};
    cell_parse_info(&cell);
    printf("nHeIII=%10.4e\n",cell.nHeIII);
    double Gamma[3], phHeating[3];
    calc_photoGH(spec, &cell, Gamma, phHeating);
    printf("eV2erg=%10.4e\n", eV2erg);
    printf("Gamma_HI=%10.4e, %10.4e, %10.4e\n",Gamma[0],Gamma[1],Gamma[2]);
    printf("Heating_HI=%10.4e, %10.4e, %10.4e\n",phHeating[0],phHeating[1],phHeating[2]);
    */
    
    return 0;
}
