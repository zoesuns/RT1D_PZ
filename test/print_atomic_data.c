#include "../base.h"
#include "../atomic_data.h"
#include "../test_func.h"
#include "fake_cells.h"
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
    char outdir[]="atomic_data_test/";
    char fpath[200];
    
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }
    strcpy(fpath, outdir);
    strcat(fpath, "eVarrayCen.dat");
    print_spec_related(eVarrayCen, 1, fpath);

    strcpy(fpath, outdir);
    strcat(fpath, "eVarrayLeft.dat");
    print_spec_related(eVarrayLeft, 1, fpath);

    strcpy(fpath, outdir);
    strcat(fpath, "eVarrayRight.dat");
    print_spec_related(eVarrayRight, 1, fpath);

    strcpy(fpath, outdir);
    strcat(fpath, "eVarrayWidth.dat");
    print_spec_related(eVarrayWidth, 1, fpath);



    
    /* print cross sections */
    double sig[NUM_ENERGY_BINS];
    for (int i=0; i<NUM_ENERGY_BINS; i++){
        sig[i]=sigma_HI_nu(eVarrayCen[i]);
    }
    strcpy(fpath, outdir);
    strcat(fpath,"sigma_HI.dat");
    print_spec_related(sig, 1, fpath);

    for (int i=0; i<NUM_ENERGY_BINS; i++){
        sig[i]=sigma_HeI_nu(eVarrayCen[i]);
    }
    strcpy(fpath, outdir);
    strcat(fpath,"sigma_HeI.dat");
    print_spec_related(sig, 1, fpath);

    for (int i=0; i<NUM_ENERGY_BINS; i++){
        sig[i]=sigma_HeII_nu(eVarrayCen[i]);
    }
    strcpy(fpath, outdir);
    strcat(fpath,"sigma_HeII.dat");
    print_spec_related(sig, 1, fpath);

    int numTpoints=71;
    double testTarray[numTpoints];
    double logTmin=1;
    double logTmax=8;
    double perbin=(logTmax-logTmin)/(numTpoints-1);

    FILE *fp;
    strcpy(fpath, outdir);
    strcat(fpath, "testTarray.dat");

    fp=fopen(fpath, "w+");

    for (int i=0; i<numTpoints; i++){
        testTarray[i]=pow(10, logTmin+perbin*i);
        fprintf(fp, "%10.6e\n", testTarray[i]);
    }
    fclose(fp);

    double eGamma[numTpoints];
    strcpy(fpath, outdir);
    strcat(fpath,"eGamma_HI.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        eGamma[i]=eGamma_HI(testTarray[i]);
        fprintf(fp, "%10.6e\n",eGamma[i]);
    }
    fclose(fp);

    strcpy(fpath, outdir);
    strcat(fpath,"eGamma_HeI.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        eGamma[i]=eGamma_HeI(testTarray[i]);
        fprintf(fp, "%10.6e\n",eGamma[i]);
    }
    fclose(fp);

    strcpy(fpath, outdir);
    strcat(fpath,"eGamma_HeII.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        eGamma[i]=eGamma_HeII(testTarray[i]);
        fprintf(fp, "%10.6e\n",eGamma[i]);
    }
    fclose(fp);


    double recomb[numTpoints];
    //init_load_recomb_data();

    strcpy(fpath, outdir);
    strcat(fpath,"recomb_HI.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        recomb[i]=recomb_rate(testTarray[i], HI);
        fprintf(fp, "%10.6e\n", recomb[i]);
    }
    fclose(fp);

    strcpy(fpath, outdir);
    strcat(fpath,"recomb_HeI.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        recomb[i]=recomb_rate(testTarray[i], HeI);
        fprintf(fp, "%10.6e\n", recomb[i]);
    }
    fclose(fp);

    strcpy(fpath, outdir);
    strcat(fpath,"recomb_HeI_r.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        recomb[i]=recomb_rate(testTarray[i], HeI_r);
        fprintf(fp, "%10.6e\n", recomb[i]);
    }
    fclose(fp);

    strcpy(fpath, outdir);
    strcat(fpath,"recomb_HeI_d.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        recomb[i]=recomb_rate(testTarray[i], HeI_d);
        fprintf(fp, "%10.6e\n", recomb[i]);
    }
    fclose(fp);

    strcpy(fpath, outdir);
    strcat(fpath,"recomb_HeII.dat");
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        recomb[i]=recomb_rate(testTarray[i], HeII);
        fprintf(fp, "%10.6e\n", recomb[i]);
    }
    fclose(fp);


    /* calc collisional ionization equilibrium */
    cellProp cell;
    cold_dense_neutral_cell(&cell,1,0.1);

    sprintf(fpath,"%sbremCoolingRateHII.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        //cell.xHI=recomb_rate(testTarray[i], HI) / (recomb_rate(testTarray[i], HI) + eGamma_HI(testTarray[i]));
        //double xHeIoxHeII=recomb_rate(testTarray[i], HeI)/eGamma_HeI(testTarray[i]);
        //double xHeIIIoxHeII=eGamma_HeII(testTarray[i])/recomb_rate(testTarray[i], HeII);
        //cell.xHeI=xHeIoxHeII/(xHeIoxHeII+1+xHeIIIoxHeII);
        //cell.xHeII=1./(xHeIoxHeII+1+xHeIIIoxHeII);
        cell.xHI=0;
        cell.nHe=0;
        cell.T=testTarray[i];
        cell_parse_info(&cell);
        fprintf(fp, "%10.6e\n", brem_cooling(cell.T,cell.ne,cell.nHII,cell.nHeII,cell.nHeIII)/cell.nH/cell.nH);
    }
    fclose(fp);

    sprintf(fpath,"%srecombCoolingRateHII.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        cell.xHI=0;
        cell.nHe=0;
        cell.T=testTarray[i];
        cell_parse_info(&cell);
        fprintf(fp, "%10.6e\n", recomb_cooling(cell.T,cell.ne,cell.nHII,cell.nHeII,cell.nHeIII)/cell.nH/cell.nH);
    }
    fclose(fp);

    sprintf(fpath,"%scolliIonizCoolingRateHII.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        //cell.xHI=recomb_rate(testTarray[i], HI) / (recomb_rate(testTarray[i], HI) + eGamma_HI(testTarray[i]));
        //double xHeIoxHeII=recomb_rate(testTarray[i], HeI)/eGamma_HeI(testTarray[i]);
        //double xHeIIIoxHeII=eGamma_HeII(testTarray[i])/recomb_rate(testTarray[i], HeII);
        //cell.xHeI=xHeIoxHeII/(xHeIoxHeII+1+xHeIIIoxHeII);
        //cell.xHeII=1./(xHeIoxHeII+1+xHeIIIoxHeII);
        cell.xHI=0.9;
        cell.nHe=0;
        cell.T=testTarray[i];
        cell_parse_info(&cell);
        fprintf(fp, "%10.6e\n", colli_ioniz_cooling(cell.T,cell.ne,cell.nHI,cell.nHeI,cell.nHeII)/cell.ne/cell.nHI);
    }
    fclose(fp);

    sprintf(fpath,"%scolliExcitCoolingRateHII.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        //cell.xHI=recomb_rate(testTarray[i], HI) / (recomb_rate(testTarray[i], HI) + eGamma_HI(testTarray[i]));
        //double xHeIoxHeII=recomb_rate(testTarray[i], HeI)/eGamma_HeI(testTarray[i]);
        //double xHeIIIoxHeII=eGamma_HeII(testTarray[i])/recomb_rate(testTarray[i], HeII);
        //cell.xHeI=xHeIoxHeII/(xHeIoxHeII+1+xHeIIIoxHeII);
        //cell.xHeII=1./(xHeIoxHeII+1+xHeIIIoxHeII);
        cell.xHI=0.9;
        cell.nHe=0;
        cell.T=testTarray[i];
        cell_parse_info(&cell);
        fprintf(fp, "%10.6e\n", colli_excit_cooling(cell.T,cell.ne,cell.nHI,cell.nHeI,cell.nHeII)/cell.ne/cell.nHI);
    }
    fclose(fp);

    sprintf(fpath,"%sinvCompCoolingz6.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        fprintf(fp, "%10.6e\n", inv_comp_cooling(testTarray[i],1,6));
    }
    fclose(fp);

    sprintf(fpath,"%sinvCompCoolingz7.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        fprintf(fp, "%10.6e\n", inv_comp_cooling(testTarray[i],1,7));
    }
    fclose(fp);


    sprintf(fpath,"%sinvCompCoolingz8.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        fprintf(fp, "%10.6e\n", inv_comp_cooling(testTarray[i],1,8));
    }
    fclose(fp);

    sprintf(fpath,"%sinvCompCoolingz9.dat",outdir);
    fp=fopen(fpath, "w+");
    for (int i=0; i<numTpoints; i++){
        fprintf(fp, "%10.6e\n", inv_comp_cooling(testTarray[i],1,9));
    }
    fclose(fp);


    /* cooling */
    /*
    cellProp denNeuCold, denIonCell;
    cell->nH=1
    cell->xHI
    parse_cell(testcell);
    

    */


    //free_loaded_recomb_data();
    return 0;
}
