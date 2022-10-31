#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<assert.h>
#include <gsl/gsl_errno.h>
#include"base.h"
#include"atomic_data.h"
#include"spec.h"
#include"calc_rates.h"
#include"qss.h"
#include"integrate.h"
#include"test_func.h"

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




int test(){
    double sa[]={0};
    double sf[]={2,5,8};
    int ind[3];
    printf("%d, %d, %d\n", ind[0], ind[1], ind[2]);
    find_index_in_sorted_array(ind,sf,3,sa,1);
    for (int i=0; i<3; i++){
        printf("index=%d\n", ind[i]);
    }
    return 0;
}



int main(){
    double zuni=6;
    
    /* initialize required quantities and functions */
    init_load_recomb_data();
    init_eVarray();

    /* a litte test:
    double tmp=recomb_rate(1.0e+01 , HI);
    printf("rate=%10.4e\n",tmp); */

    int iout_prev;
    /* define variables used in solving ode*/
    double y[4];
    double bkg[6]={0};
    double err[4]={1e-2,1e-2,1e-2,1e-2};
    cellProp cell;
 
    qss_system *qs;
    qs=qss_alloc(4, rate_func, NULL);
    rateParams  inParams;

    int numCells=400; 

    int numOut=8;
    double tOutList[]={ 1e5*yr2s, 3e5*yr2s, 1e6*yr2s, 3e6*yr2s, 1e7*yr2s, 3e7*yr2s, 6e7*yr2s, 1e8*yr2s};
    //double tOutList[]={ 1e8*yr2s};
    double dOutList[numCells];
    double yOutList[numOut][numCells][4];
/*
    double **yOutList;
    yOutList=(double **)cart_alloc_worker(numCells*sizeof(double*),__FILE__,__LINE__);
    for (int i=0; i<numCells; i++){
        yOutList[i]=cart_alloc(double, 4);
    }   
*/


    int *tIndice;

    /* initialize spectrum list, initially just one constant spectrum */

    outArrayBlock *specBlockIn;
    specBlockIn=output_alloc(NUM_ENERGY_BINS, 1);
    double *spec;
    spec=init_spec();
    new_powerlaw_spec(spec, 1e57, 1.5);
    specBlockIn->current=1;
    specBlockIn->t[0]=0;
    specBlockIn->y[0]=spec;

    outArrayBlock *specBlockOut;
    specBlockOut=output_alloc(NUM_ENERGY_BINS, 100);
    specBlockOut->current=1;
    specBlockOut->t[0]=0;

    
    outArrayBlock *tmpOut=output_alloc(4,100);

    printf("entering the loop\n");
    for (int k=0; k<numCells; k++){
        printf("cell #%d\n", k);
        double dist=0.01*k+0.1;
        fake_cell(&cell, dist, 0.01); 
        cell_parse_info(&cell);
        //bkg[6]={0,0,0,0,0,0}; // replace later by true cell and bkg

        copy_cell(&(inParams.cell_pass), cell);

        y[0]=cell.xHI; y[1]=cell.xHeI; y[2]=cell.xHeII; y[3]=cell.T;

        cell_parse_info(&cell);

        // the first output spec (at t=0)
        specBlockOut->current=1;
        specBlockOut->t[0]=0;
        calc_transmitted_spec(specBlockOut->y[0], specBlockIn->y[0], y, &cell);
        
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

                    update_inParams(&inParams, specBlockIn->y[0], y, bkg, zuni);
                    qss_solve_save(qs, ts, tOutList[i], y, err, &inParams, tmpOut);
                    dOutList[k]=dist;
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
                    for (it=0; it<tIndice[0]-1; it++){
                        iout_prev=tmpOut->current;
                        update_inParams(&inParams,specBlockIn->y[it], y, bkg, zuni);
                        qss_solve_save(qs, specBlockIn->t[it], specBlockIn->t[it+1], y, err, &inParams, tmpOut);

                        assert(iout_prev<specBlockOut->size);
                        save_all_spec(specBlockOut, iout_prev, inParams.absNdot_nu, tmpOut, &cell);
                    }
                    iout_prev=tmpOut->current;

                    update_inParams(&inParams,specBlockIn->y[it], y, bkg, zuni);
                    qss_solve_save(qs, specBlockIn->t[it], tOutList[0], y, err, &inParams, tmpOut);
                    
                    dOutList[k]=dist;
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

                    update_inParams(&inParams,specBlockIn->y[tIndice[i-1]], y, bkg, zuni);
                    qss_solve_save(qs, ts, tOutList[i], y, err, &inParams, tmpOut);
                    dOutList[k]=dist;
                    yOutList[i][k][0]=y[0];
                    yOutList[i][k][1]=y[1];
                    yOutList[i][k][2]=y[2];
                    yOutList[i][k][3]=y[3]; 

                    assert(iout_prev<specBlockOut->size);
                    save_all_spec(specBlockOut, iout_prev, specBlockIn->y[tIndice[i-1]], tmpOut, &cell);
                } else { // multiple spectra in the time interval
                    iout_prev=tmpOut->current;

                    update_inParams(&inParams,specBlockIn->y[tIndice[i-1]], y, bkg, zuni);
                    qss_solve_save(qs, ts, specBlockIn->t[tIndice[i-1]+1], y, err, &inParams, tmpOut);
                    dOutList[k]=dist;
                    yOutList[i][k][0]=y[0];
                    yOutList[i][k][1]=y[1];
                    yOutList[i][k][2]=y[2];
                    yOutList[i][k][3]=y[3]; 

                    assert(iout_prev<specBlockOut->size);
                    save_all_spec(specBlockOut, iout_prev, specBlockIn->y[tIndice[i-1]], tmpOut, &cell);

                    int it=0;
                    for (it=tIndice[i-1]+1; it<tIndice[i]; it++){
                        iout_prev=tmpOut->current;

                        update_inParams(&inParams,specBlockIn->y[it], y, bkg, zuni);
                        qss_solve_save(qs, specBlockIn->t[it], specBlockIn->t[it+1], y, err, &inParams, tmpOut);

                        assert(iout_prev<specBlockOut->size);
                        save_all_spec(specBlockOut, iout_prev, inParams.absNdot_nu, tmpOut, &cell);
                    }
                    iout_prev=tmpOut->current;

                    update_inParams(&inParams,specBlockIn->y[it], y, bkg, zuni);
                    qss_solve_save(qs, specBlockIn->t[it], tOutList[i], y, err, &inParams, tmpOut);
                    dOutList[k]=dist;
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


        char fname[100]="";
        if (k%(int)(numCells/10)==0){
            sprintf(fname, "%d", k);
            strcat(fname,"cell.out");
            write2file_ty(tmpOut->t, tmpOut->y, 4, tmpOut->current, fname);

            sprintf(fname, "%d", k);
            strcat(fname,".inspec");
            write2file_ty(specBlockIn->t, specBlockIn->y, NUM_ENERGY_BINS, specBlockIn->current, fname);

            sprintf(fname, "%d", k);
            strcat(fname,".outspec");
            write2file_ty(specBlockOut->t, specBlockOut->y, NUM_ENERGY_BINS, specBlockOut->current, fname);
        }

        printf("specBlockIn->current=%zu\n",specBlockIn->current);
        printf("AFTER EVOLUTION:\n");
        printf("tmpOut->current=%zu\n",tmpOut->current);
        printf("specBlockOut->current=%zu\n",specBlockOut->current);
        int nspec=tmpOut->current;
        free_block(specBlockIn);
        specBlockIn = specBlockOut;
        specBlockOut=output_alloc(NUM_ENERGY_BINS, nspec);

        free_block(tmpOut);
        tmpOut=output_alloc(4,nspec);
        printf("specBlockIn->current=%zu\n",specBlockIn->current);
    }

    for (int i=0; i<numOut; i++){
        char fname[100]="";
        sprintf(fname, "%4.2e", tOutList[i]/yr2s);
        strcat(fname,"yr.txt");
        write_1_4(dOutList, yOutList[i], numCells, fname);
    }

    qss_free(qs);
    free_loaded_recomb_data();

    return 0;
}
