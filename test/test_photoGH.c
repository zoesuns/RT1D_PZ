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

struct stat st = {0};




int main(){
    char outdir[]="photoGH_test/";

    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }


    init_eVarray();    
    init_load_recomb_data();

    double *spec;
    spec=init_spec();

    cellProp cell;

    void (*cellFunc[4]) (cellProp *cell, double dist_pMpc, double dr_pMpc);
    cellFunc[0]=cold_dense_neutral_cell;
    cellFunc[1]=hot_diffuse_ionized_cell;
    cellFunc[2]=H_only_cell;
    cellFunc[3]=He_only_cell;

    char outputNames[][100]={"cold.out", "hot.out", "H_only.out", "He_only.out"};
    char outputSpecNames[][100]={"cold.outspec", "hot.outspec", "H_only.outspec", "He_only.outspec"};
    char outn[200], outs[200];

    double G[3], H[3];

    char inputSpecName[][200]={"spectra_test/Ndot_1e57_alpha_3.0_NHI_1e+19_NHeI_1e+16_NHeII_1e+15.dat"};
    for (int i=0; i<4; i++){
        printf("new cell\n");
        cellFunc[i](&cell, 0.1, 0.01);
        for (int j=0; j<1; j++){
            read_1D_data_file(inputSpecName[j], NUM_ENERGY_BINS, spec, '#', 0);

            strcpy(outn, outdir);
            strcpy(outs, outdir);
            strcat(outn, outputNames[i]);
            strcat(outs, outputSpecNames[i]);
            calc_photoGH(spec, &cell, G, H);
            printf("%e %e %e %e %e %e\n", G[0], G[1], G[2], H[0], H[1], H[2]);
        }
    }

    free_loaded_recomb_data();
    return 0;
}
