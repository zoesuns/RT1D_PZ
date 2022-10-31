#include "../base.h"
#include "../atomic_data.h"
#include "../spec.h"
#include "../test_func.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct stat st = {0};



int main(){
    init_eVarray();
    char outdir[]="spectra_test/";
    char fpath[200];
    char fname[100];
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }
    
    double *spec=init_spec();

    double alphaList[12]={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 10, 100};
    for (int i=0; i<12; i++){
        new_powerlaw_spec(spec, 1e57, alphaList[i]);
        strcpy(fpath, outdir);
        sprintf(fname, "Ndot_1e57_alpha_%.1f.dat", alphaList[i]);
        strcat(fpath, fname);
        print_spec_related(spec, 1, fpath);
    }

    double NHIList[8]={1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22}; 
    double NHeIList[8]={1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22}; 
    double NHeIIList[8]={1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22}; 
    cellProp cell;
    cell.dr_cm=1e22;
    cell.nH=100;
    cell.nHe=8;
    double y[4];

    new_powerlaw_spec(spec, 1e57, 1.5);
    double *absSpec=init_spec();

    for (int i=0; i<8; i++){
        for (int j=0; j<8; j++){
            for (int k=0; k<8; k++){
                y[0]=NHIList[i]/cell.dr_cm/cell.nH;
                y[1]=NHeIList[j]/cell.dr_cm/cell.nHe;
                y[2]=NHeIIList[k]/cell.dr_cm/cell.nHe;
                calc_transmitted_spec(absSpec, spec, y, &cell);
                strcpy(fpath, outdir);
                sprintf(fname, "Ndot_1e57_alpha_1.5_NHI_%.0e_NHeI_%.0e_NHeII_%.0e.dat", NHIList[i], NHeIList[j], NHeIIList[k]);
                strcat(fpath, fname);
                print_spec_related(absSpec, 1, fpath);
            }
        }
    }

    new_powerlaw_spec(spec, 1e57, 3);

    for (int i=0; i<8; i++){
        for (int j=0; j<8; j++){
            for (int k=0; k<8; k++){
                y[0]=NHIList[i]/cell.dr_cm/cell.nH;
                y[1]=NHeIList[j]/cell.dr_cm/cell.nHe;
                y[2]=NHeIIList[k]/cell.dr_cm/cell.nHe;
                calc_transmitted_spec(absSpec, spec, y, &cell);
                strcpy(fpath, outdir);
                sprintf(fname, "Ndot_1e57_alpha_3.0_NHI_%.0e_NHeI_%.0e_NHeII_%.0e.dat", NHIList[i], NHeIList[j], NHeIIList[k]);
                strcat(fpath, fname);
                print_spec_related(absSpec, 1, fpath);
            }
        }
    }

    return 0;
}
