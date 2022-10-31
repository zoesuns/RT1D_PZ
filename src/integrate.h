#ifndef _INCL_INTEGRATE_
#define _INCL_INTEGRATE_

typedef struct rate_params_t{
    double *absNdot_nu;
    cellProp cell_pass;
    double *bkgRad;
    double zuni;
} rateParams;

int update_inParams(rateParams *inParams, double *spec, double *y, double *bkgRad, double zuni);

int save_all_spec(outArrayBlock *specBlockOut, int iprev, double *specIn, outArrayBlock *tmpOut, cellProp *cell);

void rate_func(double t, double *y, void *params, double *w, double *a);

#endif
