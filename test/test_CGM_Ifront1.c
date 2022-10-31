#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "../base.h"
#include "../atomic_data.h"
#include "../spec.h"
#include "../calc_rates.h"
#include "../propagate.h"
#include "../cosmology.h"
#include "fake_cells.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

struct stat st = {0};

int main(){
    init_load_recomb_data();
    init_load_cosmology_data();

    init_eVarray();

    char outdir[]="propagate_pure_H_CGM/CGM_02_04pMpc/"; //////////////
    char fpath[200];
    char fname[100];
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }

    int nc=800; //////////////// 20x10kpc(nH0.0001)+200x1kpc(nH0.01)+140x40kpc(nH0.0001)
    double bkg[nc][6];
    double zuni=7;
    int no=12;
    double tOutList[]={1e3*yr2s, 3e3*yr2s, 1e4*yr2s, 3e4*yr2s, 1e5*yr2s, 3e5*yr2s, 1e6*yr2s, 3e6*yr2s, 1e7*yr2s, 3e7*yr2s, 6e7*yr2s, 1e8*yr2s};
    double *spec;
    spec=init_spec();
    //new_powerlaw_spec(spec, 1e57, 5);
    new_zero_spec(spec);
    spec[0]=1e57/ eVarrayWidth[0];

    cellProp cell[nc];

    for (int k=0; k<20; k++){
        for (int j=0; j<6; j++){
            bkg[k][j]=0;
        }
        double dr=0.01; //////////////////
        double dist=dr*k+0.001;
        z7_cell(&cell[k], dist, dr);
        cell[k].nH=1e-4; ////////////
        cell[k].nHe=0;
        cell[k].T=1e4;
        cell_parse_info(&cell[k]);
    }

    for (int k=20; k<220; k++){
        for (int j=0; j<6; j++){
            bkg[k][j]=0;
        }
        double dr=0.001; //////////////////
        double dist=0.2+dr*(k-20)+0.001;
        z7_cell(&cell[k], dist, dr);
        cell[k].nH=1e-2; ////////////
        cell[k].nHe=0;
        cell[k].T=1e4;
        cell_parse_info(&cell[k]);
    }

    for (int k=220; k<nc; k++){
        for (int j=0; j<6; j++){
            bkg[k][j]=0;
        }
        double dr=0.01; //////////////////
        double dist=0.4+dr*(k-220)+0.001;
        z7_cell(&cell[k], dist, dr);
        cell[k].nH=1e-4; ////////////
        cell[k].nHe=0;
        cell[k].T=1e4;
        cell_parse_info(&cell[k]);
    }





    print_cell_properties(&cell[2]);
    propagate_los(nc, cell, spec, zuni, no, tOutList, bkg, outdir, "test");










    free_loaded_cosmology_data();

    free_loaded_recomb_data();


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
