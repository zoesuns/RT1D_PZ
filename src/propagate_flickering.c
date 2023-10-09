#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include "base.h"
#include "atomic_data.h"
#include "spec.h"
#include "calc_rates.h"
#include "qss.h"
#include "integrate.h"
#include "test_func.h"
#include "cosmology.h"

outArrayBlock* remove_similar_spec(outArrayBlock *specBlockIn, int trimSpec){
    //printf("LETS SHRINK SPECTRA BLOCK!!\n");
    outArrayBlock* specBlockInShrink=output_alloc(NUM_ENERGY_BINS, specBlockIn->current+1);
    specBlockInShrink->t[0]=specBlockIn->t[0];
    for (int j=0; j<NUM_ENERGY_BINS; j++){
        specBlockInShrink->y[0][j]=specBlockIn->y[0][j];
    }
    specBlockInShrink->current+=1;

    if (trimSpec==0){
        for (int i=1; i<specBlockIn->current; i++){
            int inewPrev=specBlockInShrink->current-1;
            if ((specBlockIn->t[i]-specBlockInShrink->t[inewPrev])<0.001*yr2s){ // set time resolution
                //printf("specBlockIn->t[%d]=%6.4e, specBlockInShrink->t[%d]=%6.4e\n", i, specBlockIn->t[i]/yr2s, inewPrev,specBlockInShrink->t[specBlockInShrink->current-1]/yr2s);
            } else{
                specBlockInShrink->t[specBlockInShrink->current]=specBlockIn->t[i];
                for (int j=0; j<NUM_ENERGY_BINS; j++){
                    specBlockInShrink->y[specBlockInShrink->current][j]=specBlockIn->y[i][j];
                }
                specBlockInShrink->current+=1;
            }
        }
    }
    else if (trimSpec==1){
        double tooDimHI=1e40;
        double tooDimHeII=1e40;
        for (int i=1; i<specBlockIn->current; i++){
            int inewPrev=specBlockInShrink->current-1;
            if ((specBlockIn->t[i]-specBlockInShrink->t[inewPrev])<0.001*yr2s){ // set time resolution
            } else{
                if ((fabs((specBlockIn->y[i][0]+tooDimHI)/(specBlockInShrink->y[inewPrev][0]+tooDimHI)-1)<1e-3) &&    //13.97eV
                        (fabs((specBlockIn->y[i][4]+tooDimHI)/(specBlockInShrink->y[inewPrev][4]+tooDimHI)-1)<1e-2) &&     //17.32eV
                            (fabs((specBlockIn->y[i][11]+tooDimHI)/(specBlockInShrink->y[inewPrev][11]+tooDimHI)-1)<1e-1) &&      //25.23eV
                                (fabs((specBlockIn->y[i][15]+tooDimHI)/(specBlockInShrink->y[inewPrev][15]+tooDimHI)-1)<1e-1) &&   //31.27eV
                                    (fabs((specBlockIn->y[i][26]+tooDimHeII)/(specBlockInShrink->y[inewPrev][26]+tooDimHeII)-1)<1e-3) &&      //56.47eV
                                        (fabs((specBlockIn->y[i][35]+tooDimHeII)/(specBlockInShrink->y[inewPrev][35]+tooDimHeII)-1)<1e-2) &&    //91.57586  eV
                                            (fabs((specBlockIn->y[i][50]+tooDimHeII)/(specBlockInShrink->y[inewPrev][50]+tooDimHeII)-1)<1e-1)){    //204.9946 eV
                } else {
                    specBlockInShrink->t[specBlockInShrink->current]=specBlockIn->t[i];
                    for (int j=0; j<NUM_ENERGY_BINS; j++){
                        specBlockInShrink->y[specBlockInShrink->current][j]=specBlockIn->y[i][j];
                    }
                    specBlockInShrink->current+=1;
                }
            }
        }
    }

    else {
        printf("Please choose if trim (1) spectra or not (0). Exiting... \n");
        exit(-1);
    }

    free_block(specBlockIn);
    return specBlockInShrink;
}


int find_index_in_sorted_array(int *indice, double *searchFor, int ns, double *sortedArray, int n){
    // sortedArray in ascending order!!!
    // searchFor also need to be in ascending order!!!
    printf("searchFor must be positive (non-zero) !!!\n");
    //printf("size=%zu,%zu, %zu\n",sizeof(sortedArray),sizeof(double*), sizeof(double));
    printf("n=%d\n",n);

    assert(searchFor[0]>0);
    assert(sortedArray[0]==0);
    int tmp=0;
    int i=0, j=0;
    for (i=0; i<ns; i++){
        indice[i]=tmp;
        //printf("tmp=%d\n",tmp);
        for (j=tmp; j<n; j++){
            if (sortedArray[j]>=searchFor[i]){
                indice[i]= j-1;
                tmp=j-1;
                break;
            }
            indice[i]=j;
            tmp=j;
        }
    }
    return 0;
}







int propagate_lc_los(int numCells, cellProp cellArr[], outArrayBlock *specBlockIn, double zuni, int numOut, double tOutList[], double bkg[][6],
                    char outputDir[], char prefix[], int trimSpec, int outCellGap){
    
    /* initialize required quantities and functions */

    int iout_prev;
    /* define variables used in solving ode*/
    double y[4];
    double err[4]={1e-2,1e-2,1e-2,1e-2};
    cellProp cell;
 
    qss_system *qs;
    qs=qss_alloc(4, rate_func, NULL);
    rateParams  inParams;
    double dOutList[numCells];
    double yOutList[numOut][numCells][4];


    int *tIndice;

    outArrayBlock *specBlockOut;
    specBlockOut=output_alloc(NUM_ENERGY_BINS, 100);
    specBlockOut->current=1;
    specBlockOut->t[0]=0;

    
    outArrayBlock *tmpOut=output_alloc(4,100);

    //printf("entering the loop\n");
    for (int k=0; k<numCells; k++){
        printf("numCells=%d\n",numCells);
        printf("cell #%d\n", k);

        copy_cell(&cell, cellArr[k]);
        copy_cell(&(inParams.cell_pass), cell);

        y[0]=cell.xHI; y[1]=cell.xHeI; y[2]=cell.xHeII; y[3]=cell.T;

        cell_parse_info(&cell);

        // assign the first output spec (at t=0)
        specBlockOut->current=1;
        specBlockOut->t[0]=0;
        calc_transmitted_spec(specBlockOut->y[0], specBlockIn->y[0], y, &cell);
        
        //for (int iii=0; iii<specBlockIn->current; iii++){
        //    printf("specBlockIn->t[%d]=%6.4e\n", iii, specBlockIn->t[iii]/yr2s);
        //}


        tIndice=cart_alloc(int, numOut);

        find_index_in_sorted_array(tIndice, tOutList, numOut, specBlockIn->t, specBlockIn->current);

        for (int kkk=0; kkk<numOut; kkk++){
            printf("tIndice=%d\n", tIndice[kkk]);
            printf("tspec=%e [yr]\n", specBlockIn->t[tIndice[kkk]]/yr2s);
        }

        tmpOut=output_alloc(4,100);
        //printf("calculating cell\n");

        double ts, te;
        for(int i=0; i<numOut; i++){
            // assign interval
            if (i==0){
                ts=0; //=tspec[0]
            } else {
                ts=tOutList[i-1];
            }
            te=tOutList[i];
            if (i==0){ // first time interval
                if (tIndice[i]==0){
                    iout_prev=tmpOut->current;
                    //zuni=cosmo_z_age(cosmicTime0+ts/(1e8*yr2s));
                    update_inParams(&inParams, specBlockIn->y[0], y, bkg[k], zuni);
                    qss_solve_save(qs, ts, tOutList[i], y, err, &inParams, tmpOut);
                    dOutList[k]=cell.dist_pMpc;
                    yOutList[i][k][0]=y[0];
                    yOutList[i][k][1]=y[1];
                    yOutList[i][k][2]=y[2];
                    yOutList[i][k][3]=y[3]; 

                    assert(iout_prev<specBlockOut->size);
                    save_all_spec( specBlockOut, iout_prev, inParams.absNdot_nu, tmpOut, &cell);
                    //print_block(specBlockOut);
                } 
                else { //(tIndice[0]>0)
                    int it=0;
                    for (it=0; it<tIndice[0]; it++){

                        iout_prev=tmpOut->current;
                        //update the first transmitted spec
                        if (it>0){
                            calc_transmitted_spec(specBlockOut->y[iout_prev], specBlockIn->y[it], y, &cell);
                        }

                        
                        //zuni=cosmo_z_age(cosmicTime0+ts/(1e8*yr2s));
                        update_inParams(&inParams,specBlockIn->y[it], y, bkg[k], zuni);
                        qss_solve_save(qs, specBlockIn->t[it], specBlockIn->t[it+1], y, err, &inParams, tmpOut);
                        //printf("t_end=%6.4e, tmpOut->t_end=%6.4e\n",specBlockIn->t[it+1]/yr2s,tmpOut->t[tmpOut->current-1]/yr2s);
                        assert(iout_prev<specBlockOut->size);
                        save_all_spec(specBlockOut, iout_prev, inParams.absNdot_nu, tmpOut, &cell);
                    }

                    iout_prev=tmpOut->current;
                    //update the first transmitted spec
                    if (it>0){
                        calc_transmitted_spec(specBlockOut->y[iout_prev], specBlockIn->y[it], y, &cell);
                    }

                    iout_prev=tmpOut->current;
                    //zuni=cosmo_z_age(cosmicTime0+ts/(1e8*yr2s));
                    update_inParams(&inParams,specBlockIn->y[it], y, bkg[k], zuni);
                    qss_solve_save(qs, specBlockIn->t[it], tOutList[0], y, err, &inParams, tmpOut);
                    //printf("t_end=%6.4e, tmpOut->t_end=%6.4e\n",tOutList[0]/yr2s,tmpOut->t[tmpOut->current-1]/yr2s);
                    //printf("tmpOut->t_end=%6.4e\n",tmpOut->t[tmpOut->current-2]/yr2s);
                    
                    dOutList[k]=cell.dist_pMpc;
                    yOutList[i][k][0]=y[0];
                    yOutList[i][k][1]=y[1];
                    yOutList[i][k][2]=y[2];
                    yOutList[i][k][3]=y[3]; 

                    assert(iout_prev<specBlockOut->size);
                    save_all_spec(specBlockOut, iout_prev, inParams.absNdot_nu, tmpOut, &cell);
                    //printf("output!!!\n");
                }
            }
 
            else { // i>0 later time interval
                //printf("doing t=%e -- %e [yr]\n", tOutList[i-1]/yr2s, tOutList[i]/yr2s);
                if ((tIndice[i]-tIndice[i-1])==0){

                    iout_prev=tmpOut->current;
                    //update the first transmitted spec
                    calc_transmitted_spec(specBlockOut->y[iout_prev], specBlockIn->y[tIndice[i-1]], y, &cell);

                    //zuni=cosmo_z_age(cosmicTime0+ts/(1e8*yr2s));
                    update_inParams(&inParams,specBlockIn->y[tIndice[i-1]], y, bkg[k], zuni);
                    qss_solve_save(qs, ts, tOutList[i], y, err, &inParams, tmpOut);
                    dOutList[k]=cell.dist_pMpc;
                    yOutList[i][k][0]=y[0];
                    yOutList[i][k][1]=y[1];
                    yOutList[i][k][2]=y[2];
                    yOutList[i][k][3]=y[3]; 

                    assert(iout_prev<specBlockOut->size);
                    save_all_spec(specBlockOut, iout_prev, specBlockIn->y[tIndice[i-1]], tmpOut, &cell);
                } else { // multiple spectra in the time interval

                    iout_prev=tmpOut->current;
                    //update the first transmitted spec
                    
                    calc_transmitted_spec(specBlockOut->y[iout_prev], specBlockIn->y[tIndice[i-1]], y, &cell);
                    

                    //zuni=cosmo_z_age(cosmicTime0+ts/(1e8*yr2s));
                    update_inParams(&inParams,specBlockIn->y[tIndice[i-1]], y, bkg[k], zuni);
                    qss_solve_save(qs, ts, specBlockIn->t[tIndice[i-1]+1], y, err, &inParams, tmpOut);
                    dOutList[k]=cell.dist_pMpc;
                    yOutList[i][k][0]=y[0];
                    yOutList[i][k][1]=y[1];
                    yOutList[i][k][2]=y[2];
                    yOutList[i][k][3]=y[3]; 

                    assert(iout_prev<specBlockOut->size);
                    save_all_spec(specBlockOut, iout_prev, specBlockIn->y[tIndice[i-1]], tmpOut, &cell);

                    int it=0;
                    for (it=tIndice[i-1]+1; it<tIndice[i]; it++){

                        iout_prev=tmpOut->current;
                        //update the first transmitted spec
                        if (it>0){
                            calc_transmitted_spec(specBlockOut->y[iout_prev], specBlockIn->y[it], y, &cell);
                        }

                        //zuni=cosmo_z_age(cosmicTime0+ts/(1e8*yr2s));
                        update_inParams(&inParams,specBlockIn->y[it], y, bkg[k], zuni);
                        qss_solve_save(qs, specBlockIn->t[it], specBlockIn->t[it+1], y, err, &inParams, tmpOut);

                        assert(iout_prev<specBlockOut->size);
                        save_all_spec(specBlockOut, iout_prev, inParams.absNdot_nu, tmpOut, &cell);
                    }

                    iout_prev=tmpOut->current;
                    //update the first transmitted spec
                    if (it>0){
                        calc_transmitted_spec(specBlockOut->y[iout_prev], specBlockIn->y[it], y, &cell);
                    }
                    //zuni=cosmo_z_age(cosmicTime0+ts/(1e8*yr2s));
                    update_inParams(&inParams,specBlockIn->y[it], y, bkg[k], zuni);
                    qss_solve_save(qs, specBlockIn->t[it], tOutList[i], y, err, &inParams, tmpOut);
                    dOutList[k]=cell.dist_pMpc;
                    yOutList[i][k][0]=y[0];
                    yOutList[i][k][1]=y[1];
                    yOutList[i][k][2]=y[2];
                    yOutList[i][k][3]=y[3]; 

                    assert(iout_prev<specBlockOut->size);
                    save_all_spec(specBlockOut, iout_prev, inParams.absNdot_nu, tmpOut, &cell);
                    //printf("output!!!\n");
                }
            } 
        }


        char fname[200];

        if (k%outCellGap==0){
            sprintf(fname, "%s/%s_%d_cell.out", outputDir, prefix, k);
            printf("%s\n",fname);
            write2file_ty(tmpOut->t, tmpOut->y, 4, tmpOut->current, fname);
           
        }

        printf("specBlockIn->current=%zu\n",specBlockIn->current);
        printf("AFTER EVOLUTION:\n");
        printf("tmpOut->current=%zu\n",tmpOut->current);
        printf("specBlockOut->current=%zu\n",specBlockOut->current);
        int nspec=tmpOut->current;
        // remove prevoius incidental spectra
        free_block(specBlockIn);

        // assign new incidental spectra to the next cell
        specBlockIn = remove_similar_spec(specBlockOut, trimSpec);
        // allocate new transmitted spectra block
        specBlockOut=output_alloc(NUM_ENERGY_BINS, nspec);
        free_block(tmpOut);
        tmpOut=output_alloc(4,nspec);
        printf("specBlockIn->current=%zu\n",specBlockIn->current);
    }

    for (int i=0; i<numOut; i++){
        char fname[100]="";
        sprintf(fname, "%s/%s_%4.2eyr.txt", outputDir, prefix, tOutList[i]/yr2s);
        write_1_4(dOutList, yOutList[i], numCells, fname);
    }

    qss_free(qs);

    return 0;
}
