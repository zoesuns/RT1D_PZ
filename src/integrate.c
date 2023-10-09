#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include "base.h"
#include "atomic_data.h"
#include "spec.h"
#include "calc_rates.h"
#include "qss.h"
#include "integrate.h"
#include "cosmology.h"

int update_inParams(rateParams *inParams, double *spec, double *y, double *bkgRad, double zuni){
    inParams->cell_pass.xHI=y[0];
    inParams->cell_pass.xHeI=y[1];
    inParams->cell_pass.xHeII=y[2];
    inParams->cell_pass.T=y[3];
    cell_parse_info(&(inParams->cell_pass));
    inParams->absNdot_nu=spec;
    inParams->bkgRad=bkgRad;
    inParams->zuni=zuni;
    return 0;
}

int save_all_spec(outArrayBlock *specBlockOut, int iprev, double *specIn, outArrayBlock *tmpOut, cellProp *cell){
    //printf("inside save_all_spec, iprev=%d, tmpOut->current=%d\n", iprev, tmpOut->current);
    for(int jj=iprev; jj<tmpOut->current; jj++){
        if (specBlockOut->size<(1+tmpOut->current)){
            output_alloc_extend(specBlockOut, NUM_ENERGY_BINS);
        }
        specBlockOut->t[jj+1]=tmpOut->t[jj];
        calc_transmitted_spec(specBlockOut->y[jj+1], specIn, tmpOut->y[jj], cell);
        specBlockOut->current+=1;
    }
    return 0;
}


void rate_func(double t, double *y, void *params, double *w, double *a){
    // dy/dt = w-a*y
    // params should contain info about
    // spectra absNdot_nu
    // distance, dr

    if (y[3]<10.) {y[3]=10;} //set a T floor here

    cellProp cell_init = ((rateParams *) params)->cell_pass;

    cellProp newCell=cell_init;
    newCell.xHI=y[0]; newCell.xHeI=y[1]; newCell.xHeII=y[2]; newCell.T=y[3];
    cell_parse_info(&newCell);

    double Gamma_qso[3], phHeating_qso[3], Cooling;
    double cosmicTime0=cosmo_age_z(((rateParams *) params)->zuni);
    double zuni_now=cosmo_z_age(cosmicTime0+t/(1e6*yr2s));
    //printf("cosmicTime0=%f,zuni_now=%f\n",cosmicTime0,zuni_now);
    //print_cell_properties(&newCell);
    Cooling=calc_Cooling(&newCell, zuni_now);
    //printf("check point 0\n");
    double Hubble_cgs=cosmo_HubbleConst_z(zuni_now)*1e5/Mpc2cm; // need to change later
    //printf("t=%e, zuni_0=%e, zuni_now=%e", t, ((rateParams *) params)->zuni, zuni_now);

    double Gamma_bkg[3]={((rateParams *) params)->bkgRad[0], ((rateParams *) params)->bkgRad[1], ((rateParams *) params)->bkgRad[2]};
    double phHeating_bkg[3]={((rateParams *) params)->bkgRad[3], ((rateParams *) params)->bkgRad[4], ((rateParams *) params)->bkgRad[5]};

    calc_photoGH(((rateParams *) params)->absNdot_nu, &newCell, Gamma_qso, phHeating_qso);
    //printf("check point 1\n");

    double alpha[3]={recomb_rate(newCell.T, HI), recomb_rate(newCell.T, HeI), recomb_rate(newCell.T, HeII)};
    double eGamma[3]={eGamma_HI(newCell.T), eGamma_HeI(newCell.T), eGamma_HeII(newCell.T)};


    w[0]= newCell.xHII*newCell.ne*alpha[0];
    a[0]= Gamma_qso[0] +Gamma_bkg[0] +newCell.ne*eGamma[0];
    w[1]= newCell.xHeII*newCell.ne*alpha[1];
    a[1]= Gamma_qso[1] +Gamma_bkg[1] +newCell.ne*eGamma[1];

    w[2]= newCell.xHeI*a[1]+ alpha[2]*(1-newCell.xHeI)*newCell.ne;
    a[2]= Gamma_qso[2] +Gamma_bkg[2] +newCell.ne*eGamma[2]
            +newCell.ne*alpha[2]+newCell.ne*alpha[1];
 
    double dntotdt = -newCell.nH*(w[0]-a[0]*y[0]) -2*newCell.nHe*(w[1]-a[1]*y[1]) -newCell.nHe*(w[2]-a[2]*y[2]);
    w[3]= (2./(3.*k_Boltzmann*newCell.ntot))*
                    (newCell.nHI*(phHeating_qso[0]+phHeating_bkg[0]) 
                     +newCell.nHeI*(phHeating_qso[1]+phHeating_bkg[1])
                     +newCell.nHeII*(phHeating_qso[2]+phHeating_bkg[2])-Cooling); 
    a[3]= 2.*Hubble_cgs+dntotdt/newCell.ntot;


//    if ((cell_init.dist_pMpc>4.36807 )&&(t>1e7*yr2s)&&(cell_init.dist_pMpc<4.36910)){
//        printf("dist=%.7f pMpc, t=%.7e yr\n", cell_init.dist_pMpc, t/yr2s);
//        printf("%e %e %e %e %e %e %e %e\n",a[0], a[1], a[2], a[3],w[0], w[1], w[2], w[3]);
        //printf("a=%e %e %e %e\n",a[0], a[1], a[2], a[3]);
        //printf("w=%e %e %e %e\n",w[0], w[1], w[2], w[3]);
//    }
    if (w[0]!=w[0] || w[1]!=w[1] || w[2]!=w[2] || w[3]!=w[3] || a[0]!=a[0] || a[1]!=a[1] || a[2]!=a[2] || a[3]!=a[3]) {
        printf("newCell.xHI=%e, newCell.xHeI=%e, newCell.xHeII=%e, newCell.T=%e\n", newCell.xHI, newCell.xHeI, newCell.xHeII, newCell.T);
        printf("NAN APPEARED!!!\n");
        printf("NAN APPEARED!!!\n");
        printf("NAN APPEARED!!!\n");
        exit(-1);
    }
}

int test_integrate_one_cell(){
    init_load_recomb_data();
    double tmp=recomb_rate(1.0e+01 , HI);
    printf("rate=%10.4e\n",tmp);

 
    qss_system *qs;
    
    double err[4]={1e-1,1e-1,1e-1,1e-1};
    //double err[4]={1e-2,1e-2,1e-2,1e-2};

    cellProp cell={1, 0.01, 0,0,
                    1e-4,1e-5,
                    1,0,0.,1e3,
                    0,0, 0,0, 0,0,0,0,0};
    cell_parse_info(&cell);

    double y[4]={cell.xHI, cell.xHeI, cell.xHeII, cell.T};

    init_eVarray();

    double *spec;
    spec=init_spec();

    new_powerlaw_spec(spec, 1e57, 1.5);
    
    double bkg[6]={0,0,0,0,0,0};
    
    rateParams  inParams={spec, cell, bkg, 6};


    qs=qss_alloc(4, rate_func, NULL);
    printf("there are %d eqn\n",qs->num_eqn);

    outArrayBlock *tmpOut=output_alloc(4,100);

    double tout[3]={0,1e7*yr2s,1e8*yr2s};
    printf("y_0=%10.4e,%10.4e,%10.4e,%10.4e\n",y[0],y[1],y[2],y[3]);
    for (int j=0; j<2; j++){
        printf("j=%d\n",j);
        qss_solve_save(qs,tout[j],tout[j+1],y,err, &inParams, tmpOut);
        inParams.cell_pass=cell;
        inParams.cell_pass.xHI=y[0];
        inParams.cell_pass.xHeI=y[1];
        inParams.cell_pass.xHeII=y[2];
        inParams.cell_pass.T=y[3];
        cell_parse_info(&(inParams.cell_pass));
        printf("y_1=%10.4e,%10.4e,%10.4e,%10.4e\n",y[0],y[1],y[2],y[3]);

    }

    for (int k=0; k<tmpOut->current; k++){
        printf("t=%10.4e [yr]\n", tmpOut->t[k]/yr2s);
        printf("y_1=%10.4e,%10.4e,%10.4e,%10.4e\n",tmpOut->y[k][0],tmpOut->y[k][1],
                    tmpOut->y[k][2],tmpOut->y[k][3]);
    }
    printf("WORK?");
    write2file_ty(tmpOut->t, tmpOut->y, 4, tmpOut->current-1, "output.out");
    qss_free(qs);
    free_loaded_recomb_data();
    return 0;
}
