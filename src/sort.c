#include <math.h>
#include <stdio.h>
#include <assert.h>

int find_index_in_sortedArray(int *indice, double *searchFor, int ns, double *sortedArray, int n){
    // sortedArray in ascending order!!!
    // searchFor also need to be in ascending order!!!
    printf("searchFor must be positive (non-zero) !!!\n");
    printf("size=%zu,%zu, %zu\n",sizeof(sortedArray),sizeof(double*), sizeof(double));
    printf("n=%d\n",n);

    assert(searchFor[0]>0);
    assert(sortedArray[0]==0);
    int tmp=0;
    int i=0, j=0;
    for (i=0; i<ns; i++){
        indice[i]=tmp;
        printf("tmp=%d\n",tmp);
        for (j=tmp; j<n; j++){
            if (sortedArray[j]>=searchFor[i]){
                printf("yes\n");
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

int solve(double a, double b, double c){
    printf("solve t = %e -- %e [yr] using spectrum at t=%e\n", a, b, c);
    return 0;
}

int savespec(){
    return 0;
}

int main(){
    int ns=5;
    int nspec=7;
    double tspec[]={0, 2e4, 3e4, 4e4, 5e6, 1e7, 7e7};
    double output[]={1e5, 1e6, 1e7, 3e7, 1e8};
    int ind[ns];
    printf("size ind=%zu\n",sizeof(ind));
    find_index_in_sortedArray(ind, output,ns,tspec,nspec);
    for (int i=0; i<ns; i++){
        printf("index=%d\n", ind[i]);
    }

    assert(tspec[0]==0);
    assert(output[0]>0);
    double ts, te;
    for (int i=0; i<ns; i++){
        // assign interval
        if (i==0){
            ts=0; //=tspec[0]
        } else {
            ts=output[i-1];
        }
        te=output[i];
        //
        if (i==0){
            if (ind[i]==0){
                solve (ts, output[i], tspec[0]);
            } else { //(ind[0]>0)
                int k=0;
                for (k=0; k<ind[i]; k++){
                    solve ( tspec[k], tspec[k+1], tspec[k]);
                    savespec();
                }
                solve(tspec[k], output[0], tspec[k]);
                printf("output!!!\n");
            }
        } else { // i>0
            if ((ind[i]-ind[i-1])==0){
                printf("yes?");
                printf("ts=%e\n",ts);
                solve (ts, output[i], tspec[0]);
            } else {
                int k=0;
                solve (ts, tspec[ind[i-1]+1], tspec[ind[i-1]]);
                savespec();
                for (k=ind[i-1]+1; k<ind[i]; k++){
                    solve ( tspec[k], tspec[k+1], tspec[k]);
                    savespec();
                }
                solve(tspec[k], output[i], tspec[k]);
                printf("output!!!\n");
            }
        } 
    }

    return 0;
}



