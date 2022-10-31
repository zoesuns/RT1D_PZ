#include "base.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


extern int errno ;

int copy_cell(cellProp *copy, cellProp original){
    copy->dist_pMpc = original.dist_pMpc;
    copy->dr_pMpc = original.dr_pMpc;
    copy->dist_cm = original.dist_cm;
    copy->dr_cm = original.dr_cm;
    copy->nH = original.nH;
    copy->nHe = original.nHe;
    copy->xHI = original.xHI;
    copy->xHeI = original.xHeI;
    copy->xHeII = original.xHeII;
    copy->T = original.T;
    copy->xHI = original.xHI;
    copy->xHeI = original.xHeI;
    copy->xHeII = original.xHeII;
    copy->nHI = original.nHI;
    copy->nHII = original.nHII;
    copy->nHeI = original.nHeI;
    copy->nHeII = original.nHeII;
    copy->nHeIII = original.nHeIII;
    copy->ne = original.ne;
    copy->ntot = original.ntot;
    return 0;
}

int cell_parse_info(cellProp *cell){
    cell->dist_cm=cell->dist_pMpc*Mpc2cm;
    cell->dr_cm=cell->dr_pMpc*Mpc2cm;

    if (cell->xHI<0) {cell->xHI=1e-30;}
    if (cell->xHeI<0) {cell->xHeI=1e-30;}
    if (cell->xHeII<0) {cell->xHeII=1e-30;}
    if (cell->xHI>1) {cell->xHI=1;}
    if (cell->xHeI>1) {cell->xHeI=1;}
    if (cell->xHeII>1) {cell->xHeII=1;}
    if (cell->T<10) {cell->T=10;}

    cell->xHII=1-cell->xHI;
    cell->xHeIII=1-cell->xHeI-cell->xHeII;

    cell->nHI=cell->nH*cell->xHI;
    cell->nHII=cell->nH*cell->xHII;
    cell->nHeI=cell->nHe*cell->xHeI;
    cell->nHeII=cell->nHe*cell->xHeII;
    cell->nHeIII=cell->nHe*cell->xHeIII;

    cell->ne=cell->nHII+cell->nHeII+cell->nHeIII*2;
    cell->ntot=cell->nH+cell->nHe+cell->ne;
    return 0;
}


void print_cell_properties(cellProp *cell){
    printf("cell properties\n");
    printf("dist= %10.4e pMpc, dr= %10.4e pMpc \n", cell->dist_pMpc, cell->dr_pMpc);
    printf("dist= %10.4e   cm, dr= %10.4e cm \n", cell->dist_cm, cell->dr_cm);
    printf("ionic fraction\n");
    printf(" xHI= %10.4e,  xHeI= %10.4e\n", cell->xHI, cell->xHeI);
    printf("xHII= %10.4e, xHeII= %10.4e\n", cell->xHII, cell->xHeII);
    printf("                 xHeIII= %10.4e\n", cell->xHeIII);
    printf("densities [cm^-3]\n");
    printf("  nH = %10.4e,    nHe = %10.4e\n", cell->nH, cell->nHe);
    printf(" nHI = %10.4e,   nHeI = %10.4e\n", cell->nHI, cell->nHeI);
    printf("nHII = %10.4e,  nHeII = %10.4e\n", cell->nHII, cell->nHeII);
    printf("                   nHeIII = %10.4e\n", cell->nHeIII);
    printf("ne=%10.4e\n\n",cell->ne);
    printf("T=%10.4eK\n",cell->T);
}



int write2file_ty(double *t, double **y, int num_var, int length, char *fileName){
        FILE * fp;
        fp=fopen(fileName,"w+");
        for (int i=0; i<length; i++){
            fprintf(fp, "%10.6e ", t[i]);
            for(int j=0; j<num_var-1; j++){
                fprintf(fp, "%10.6e ", y[i][j]);
            }
            fprintf(fp, "%10.6e\n", y[i][num_var-1]);
        }
        fclose(fp);
    return 0;
}



int write_1_4(double d[], double y[][4], int length, char *fileName){
    FILE * fp;
    fp=fopen(fileName,"w+");
    int i=0, j=0;
    for (i=0; i<length; i++){
        fprintf(fp, "%10.6e ", d[i]);
        for(j=0; j<4; j++){
            fprintf(fp, "%10.6e ", y[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}



//int init_eVarray(double *nu, double *nuL, double *nuR){
int init_eVarray(){
    double logmin=log10(NU_MIN);
    double logmax=log10(NU_MAX);
    double perbin=(logmax-logmin)/NUM_ENERGY_BINS;
    for (int i=0; i<NUM_ENERGY_BINS; i++){
        eVarrayLeft[i]=pow(10,logmin+((double)i)*perbin);
        eVarrayCen[i]=pow(10,logmin+((double)i+0.5)*perbin);
        eVarrayRight[i]=pow(10,logmin+((double)i+1.)*perbin);
        eVarrayWidth[i]=eVarrayRight[i]-eVarrayLeft[i];
    }
    return 0;
}

int print_spec_related(double *eV, int ifWrite, char *fileName){
    for (int i=0; i<NUM_ENERGY_BINS; i++)
        printf("%10.6e\n",eV[i]);
    if (ifWrite){
        FILE * fp;
        fp=fopen(fileName,"w+");
        for (int i=0; i<NUM_ENERGY_BINS; i++)
            fprintf(fp, "%10.6e\n",eV[i]);
        fclose(fp);
    }
    return 0;
}

int count_nrow_data_file(char *filename, char headerMark, int *nh){
    FILE* fp;
    char chr;
    int nrow=0;
    *nh=0;
    if((fp = fopen(filename, "r")) == NULL){
        printf("%s does not exist!!!\n Exiting ... \n", filename);
        fprintf(stderr, "Value of errno: %d\n", errno);
        perror("Error printed by perror");
        fprintf(stderr, "Error opening file: %s\n", strerror( errno ));
        exit(-1);
    }
    chr=getc(fp);
    while(chr!=EOF){
        //printf("chr=%c ",chr);
        //printf("headerMark=%c ", headerMark);
        if (chr==headerMark)
            (*nh)+=1;
        else if (chr=='\n')
            nrow+=1;
        chr=getc(fp);
        //printf("nrow=%d\n",nrow);
        //printf("nh=%d\n",*nh);
        
    }
    nrow-=(*nh);
    fclose(fp);
    return nrow;
}


int read_2D_data_file(char *filename, int nrow, double *x, double *y, char headerMark, int nh)
{
    FILE* fp;
    char chr,str1[100], str2[100];
    if((fp = fopen(filename, "r")) == NULL){
        printf("%s does not exist!!!\n Exiting ... \n", filename);
        fprintf(stderr, "Value of errno: %d\n", errno);
        perror("Error printed by perror");
        fprintf(stderr, "Error opening file: %s\n", strerror( errno ));
        exit(-1);
    }
    chr=getc(fp);
    int cnt=0;
    while(chr!=EOF){
        while (chr==headerMark){
            cnt+=1;
            while (chr!='\n'){
                chr=getc(fp);
            }
            chr=getc(fp);
        }
        if (cnt==nh)
            break;
    }

    chr=ungetc(chr, fp);


    for (int i=0; i<nrow; i++){
        fscanf(fp, "%s %s", str1, str2);
        x[i]=strtod(str1,NULL);
        y[i]=strtod(str2,NULL);
        chr=getc(fp);
    }
    return 0;
}

int read_1D_data_file(char *filename, int nrow, double *y, char headerMark, int nh)
{
    FILE* fp;
    char chr,str1[100];
    if((fp = fopen(filename, "r")) == NULL){
        printf("%s does not exist!!!\n Exiting ... \n", filename);
        fprintf(stderr, "Value of errno: %d\n", errno);
        perror("Error printed by perror");
        fprintf(stderr, "Error opening file: %s\n", strerror( errno ));
        exit(-1);
    }
    chr=getc(fp);
    int cnt=0;
    while(chr!=EOF){
        while (chr==headerMark){
            cnt+=1;
            while (chr!='\n'){
                chr=getc(fp);
            }
            chr=getc(fp);
        }
        if (cnt==nh)
            break;
    }

    chr=ungetc(chr, fp);

    for (int i=0; i<nrow; i++){
        fscanf(fp, "%s", str1);
        y[i]=strtod(str1,NULL);
        chr=getc(fp);
    }
    return 0;
}






int write_1D_array_to_txt(double a1[], int nrow, char *filename, char* fmt, char* header){
    FILE* fp;
    char format[100];
    strcpy(format,fmt);
    fp=fopen(filename,"w+");
    if (strcmp(header,"") != 0){
        fprintf(fp, "%s\n", header);
    }
    if (strcmp(format,"") == 0){
        strcpy(format, "%e");
    } 
    for (int i=0; i<nrow; i++){
        fprintf(fp, format, a1[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}
    
int read_1D_array_from_txt(char *filename, int nrow, double a1[], char headerMark){
    FILE* fp;
    char chr;
    if((fp = fopen(filename, "r")) == NULL){
        printf("%s does not exist!!!\n Exiting ... \n", filename);
        fprintf(stderr, "Value of errno: %d\n", errno);
        perror("Error printed by perror");
        fprintf(stderr, "Error opening file: %s\n", strerror( errno ));
        exit(-1);
    }

    chr=getc(fp);
    //printf("chr=%c\n", chr);

    while (chr==headerMark){
        //printf("chr=%c\n", chr);
        while (chr!='\n'){
            //printf("chr=%c\n", chr);
            chr=getc(fp);
        }
        chr=getc(fp);
    }

    ungetc(chr, fp);

    for (int i=0; i<nrow; i++){
        fscanf(fp, "%le\n", &a1[i]);
    }
    return 0;
}


int write_2D_array_to_txt(double a1[], double a2[], int nrow, char *filename, char* fmt, char* header){
    FILE* fp;
    char format[100];
    strcpy(format,fmt);
    fp=fopen(filename,"w+");
    if (strcmp(header,"") != 0){
        fprintf(fp, "%s\n", header);
    }
    if (strcmp(format,"") == 0){
        strcpy(format, "%e %e");
    } 
    for (int i=0; i<nrow; i++){
        fprintf(fp, format, a1[i], a2[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}

int read_2D_array_from_txt(char *filename, int nrow, double a1[], double a2[], char headerMark){
    FILE* fp;
    char chr;
    if((fp = fopen(filename, "r")) == NULL){
        printf("%s does not exist!!!\n Exiting ... \n", filename);
        fprintf(stderr, "Value of errno: %d\n", errno);
        perror("Error printed by perror");
        fprintf(stderr, "Error opening file: %s\n", strerror( errno ));
        exit(-1);
    }

    chr=getc(fp);
    //printf("chr=%c\n", chr);

    while (chr==headerMark){
        
        //printf("chr=%c\n", chr);
        while (chr!='\n'){
            //printf("chr=%c\n", chr);
            chr=getc(fp);
        }
        chr=getc(fp);
    }

    ungetc(chr, fp);

    for (int i=0; i<nrow; i++){
        fscanf(fp, "%le %le\n", &a1[i], &a2[i]);
    }
    return 0;
}


int read_array_from_txt(char *filename, int nrow, char headerMark, int n_args, ...){
        
    double array[nrow][n_args];

    FILE* fp;
    char chr;
    if((fp = fopen(filename, "r")) == NULL){
        printf("%s does not exist!!!\n Exiting ... \n", filename);
        fprintf(stderr, "Value of errno: %d\n", errno);
        perror("Error printed by perror");
        fprintf(stderr, "Error opening file: %s\n", strerror( errno ));
        exit(-1);
    }

    chr=getc(fp);
    //printf("chr=%c\n", chr);

    while (chr==headerMark){
        
        //printf("chr=%c\n", chr);
        while (chr!='\n'){
            //printf("chr=%c\n", chr);
            chr=getc(fp);
        }
        chr=getc(fp);
    }

    ungetc(chr, fp);

    for (int i=0; i<nrow; i++){
        for (int j=0; j<n_args; j++){
            fscanf(fp, "%le\n", &array[i][j]);
        }
        fscanf(fp, "\n");
    }

    register int i;
    va_list ap;
    double* a;

    va_start(ap, n_args);

    for (i=1; i<=n_args; i++){
        a = va_arg(ap, double*);
        for (int j=0; j<nrow; j++){
            a[j]=array[j][i-1];
        }
    }

    va_end(ap);
    return 0;
}


