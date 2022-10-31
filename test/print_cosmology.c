#include "../base.h"
#include "../atomic_data.h"
#include "../test_func.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../cosmology.h"
#include "../setup_path.h"

struct stat st = {0};



int main(){
    init_eVarray();
    char outdir[]="cosmology_test/";
    char fpath[200];
    
    printf("%s", dataDir);
    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }
    sprintf(fpath, "%sH.txt", outdir);

    init_load_cosmology_data();

    double zl[100];
    double agel[100];
    for (int i=0; i<100; i++){
        zl[i]=19-14*(i/100.);
    }

    for (int i=0; i<100; i++){
        agel[i]=200+900*(i/100.);
    }

    double H[100], rhoc[100], oz[100], oage[100];
    for (int i=0; i<100; i++){
        H[i]=cosmo_HubbleConst_z(zl[i]);
        rhoc[i]=cosmo_rhoc_z(zl[i]);
        oage[i]=cosmo_age_z(zl[i]);
        oz[i]=cosmo_z_age(agel[i]);
    }
    sprintf(fpath,"%sH_z.txt",outdir);
    write_2D_array_to_txt(zl,H,100,fpath,"%e %e","#z H");
    sprintf(fpath,"%srhoc_z.txt",outdir);
    write_2D_array_to_txt(zl,rhoc,100,fpath,"%e %e","#z rhoc");
    sprintf(fpath,"%sage_z.txt",outdir);
    write_2D_array_to_txt(zl,oage,100,fpath,"%e %e","#z age[Myr]");
    sprintf(fpath,"%sz_age.txt",outdir);
    write_2D_array_to_txt(agel,oz,100,fpath,"%e %e","#age[Myr] z");




    free_loaded_cosmology_data();
    return 0;
}
