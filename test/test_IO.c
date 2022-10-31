#include "../atomic_data.h"
#include "../base.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct stat st = {0};

int main(){
    char outdir[]="IO_test/";


    if (stat(outdir, &st) == -1) {
        mkdir(outdir, 0700);
    }

/*

    double testArr[10]={1,2,3,4,5,6,7,5.5,90,2};
    double out[10];
    char fpath[200];
    strcpy(fpath, outdir);
    strcat(fpath, "test1D.txt");
    printf("here?\n");
    write_1D_array_to_txt(testArr, 10, fpath, "%2.3e", "");
    read_1D_array_from_txt(fpath, 10, out, '#');
    for (int i=0; i<10; i++){
        printf("%e\n",out[i]);
    }
*/
    double a[4], b[4], c[4];
    read_array_from_txt("IO_test/ndarray.txt", 4, '#', 3, a, b, c);
    for (int i=0; i<4; i++){
        printf("%e %e %e\n", a[i], b[i], c[i]);
    }
    return 0;
}
