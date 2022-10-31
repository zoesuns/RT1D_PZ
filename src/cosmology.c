#include <math.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "setup_path.h"
#include "atomic_data.h"
#include "base.h"

#define NUM_Z_SAMPLE 10001

//extern char dataDir[];
gsl_spline *cosmo_z_age_spline;
gsl_spline *cosmo_age_z_spline;
gsl_spline *cosmo_HubbleConst_z_spline;
gsl_spline *cosmo_rhoc_z_spline;
gsl_spline *cosmo_Om_z_spline;
gsl_spline *cosmo_Ob_z_spline;
gsl_interp_accel *acc2;

int reverse(double a[], int num){
    // https://www.cs.cmu.edu/~rbd/papers/cmj-float-to-int.html
    double tmp;
    for (int i=0; i<(int)(num/2); i++){
        tmp=a[i];
        a[i]=a[num-1-i];
        a[num-1-i]=tmp;
    }
    return 0;
}

int init_load_cosmology_data(){
    acc2 = gsl_interp_accel_alloc ();
    char fname[200];
    double z[NUM_Z_SAMPLE], age[NUM_Z_SAMPLE], HubbleConst[NUM_Z_SAMPLE], rhoc[NUM_Z_SAMPLE];
    
    sprintf(fname, "%scosmology_age_z.dat", dataDir);
    read_2D_array_from_txt(fname, NUM_Z_SAMPLE, z, age, '#');
    cosmo_z_age_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_Z_SAMPLE);
    gsl_spline_init (cosmo_z_age_spline, age, z, NUM_Z_SAMPLE);
    printf("done cosmo_z_age_spline\n");
    cosmo_age_z_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_Z_SAMPLE);
    reverse(z,NUM_Z_SAMPLE); reverse(age,NUM_Z_SAMPLE);
    for (int i=0; i<NUM_Z_SAMPLE-1; i++){
        if(!(z[i]<z[i+1])){
        printf("%e %e\n", z[i], z[i+1]);}
    }
    gsl_spline_init (cosmo_age_z_spline, z, age, NUM_Z_SAMPLE);
    printf("done cosmo_age_z_spline\n");



    sprintf(fname, "%scosmology_Hz_z.dat",dataDir);
    read_2D_array_from_txt(fname, NUM_Z_SAMPLE, z, HubbleConst, '#');
    cosmo_HubbleConst_z_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_Z_SAMPLE);
    reverse(z,NUM_Z_SAMPLE); reverse(HubbleConst,NUM_Z_SAMPLE);
    gsl_spline_init (cosmo_HubbleConst_z_spline, z, HubbleConst, NUM_Z_SAMPLE);
    printf("done cosmo_HubbleConst_z_spline\n");

    sprintf(fname, "%scosmology_rhoc_z.dat",dataDir);
    read_2D_array_from_txt(fname, NUM_Z_SAMPLE, z, rhoc, '#');
    cosmo_rhoc_z_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_Z_SAMPLE);
    reverse(z,NUM_Z_SAMPLE); reverse(rhoc,NUM_Z_SAMPLE); 
    gsl_spline_init (cosmo_rhoc_z_spline, z, rhoc, NUM_Z_SAMPLE);

    double Om[NUM_Z_SAMPLE], Ob[NUM_Z_SAMPLE];
    sprintf(fname, "%scosmology_Om_z.dat",dataDir);
    read_2D_array_from_txt(fname, NUM_Z_SAMPLE, z, Om, '#');
    cosmo_Om_z_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_Z_SAMPLE);
    reverse(z,NUM_Z_SAMPLE); reverse(Om,NUM_Z_SAMPLE); 
    gsl_spline_init (cosmo_Om_z_spline, z, Om, NUM_Z_SAMPLE);

    sprintf(fname, "%scosmology_Ob_z.dat",dataDir);
    read_2D_array_from_txt(fname, NUM_Z_SAMPLE, z, Ob, '#');
    cosmo_Ob_z_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_Z_SAMPLE);
    reverse(z,NUM_Z_SAMPLE); reverse(Ob,NUM_Z_SAMPLE); 
    gsl_spline_init (cosmo_Ob_z_spline, z, Ob, NUM_Z_SAMPLE);

    return 0;
}

double cosmo_age_z(double z){
    assert(z<19.9 && z>4.9961);
    return gsl_spline_eval (cosmo_age_z_spline, z, acc2);
}

double cosmo_z_age(double age){
    assert(age>180.5 && age<1180.4);
    return gsl_spline_eval (cosmo_z_age_spline, age, acc2);
}

double cosmo_HubbleConst_z(double z){
    assert(z<19.9 && z>4.9961);
    return gsl_spline_eval (cosmo_HubbleConst_z_spline, z, acc2);
}

double cosmo_rhoc_z(double z){
    assert(z<19.9 && z>4.9961);
    return gsl_spline_eval (cosmo_rhoc_z_spline, z, acc2);
}

double cosmo_Om_z(double z){
    assert(z<19.9 && z>4.9961);
    return gsl_spline_eval (cosmo_Om_z_spline, z, acc2);
}

double cosmo_Ob_z(double z){
    assert(z<19.9 && z>4.9961);
    return gsl_spline_eval (cosmo_Ob_z_spline, z, acc2);
}

int free_loaded_cosmology_data(){
    gsl_spline_free (cosmo_age_z_spline);
    gsl_spline_free (cosmo_z_age_spline);
    gsl_spline_free (cosmo_HubbleConst_z_spline);
    gsl_spline_free (cosmo_rhoc_z_spline);
    gsl_interp_accel_free (acc2);
    return 0;
}


