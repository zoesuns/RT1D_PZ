#include<math.h>
#include<stdio.h>
#include <gsl/gsl_errno.h>
#include"test_func.h"
#include"base.h"
#include"atomic_data.h"
#include"spec.h"
#include"calc_rates.h"
#include"qss.h"
#include"integrate.h"


void test_rate_func(double t, double *y, void *params, double *w, double *a){
    w[0]=0; 
    w[1]=2;
    w[2]=2;
    w[3]=1;
    a[0]=0;
    a[1]=0;
    a[2]=0;
    a[3]=((double *)params)[2];
}

void adjust_func(double t, double *y, void *params){
    
}


int main(){
    init_load_recomb_data();
    double tmp=recomb_rate(1.0e+01 , HI);
    printf("rate=%10.4e\n",tmp);

 
    qss_system *qs;
    
    double err[4]={1e-2,1e-2,1e-2,1e-2};
    //double err[4]={1e-2,1e-2,1e-2,1e-2};

    //cellProp cell={1, 0.01, 0,0,
    //                1e-4,1e-4,
    //                1,1,0.,1e4,
    //                0,0, 0,0, 0,0,0,0,0};
    cellProp cell;
    fake_cell(&cell, 0.1, 0.01); 

    cell_parse_info(&cell);

    double y[4]={cell.xHI, cell.xHeI, cell.xHeII, cell.T};

    init_eVarray();

    double *spec;
    spec=init_spec();

    new_powerlaw_spec(spec, 1e57, 100);
    
    double bkg[6]={0,0,0,0,0,0};
    
    rateParams  inParams={spec, cell, bkg, 6};


    qs=qss_alloc(4, rate_func, adjust_func);
    printf("there are %d eqn\n",qs->num_eqn);

    outArrayBlock *tmpOut=output_alloc(4,100);

    double tout[1]={1e8*yr2s};
    printf("y_0=%10.4e,%10.4e,%10.4e,%10.4e\n",y[0],y[1],y[2],y[3]);
    for (int j=0; j<1; j++){
        printf("j=%d\n",j);
        qss_solve_save(qs,0,tout[j],y,err, &inParams, tmpOut);
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
    write2file_ty(tmpOut->t, tmpOut->y, 4, tmpOut->current, "output.out");
    qss_free(qs);
    free_loaded_recomb_data();
    return 0;
}
