# Compiler
CC=gcc
# Compiler options, including optimization
CFLAGS= -g -Wall -std=c99
LIB= -L/data/hqchen/software/clibs/lib -lm -lgsl -lgslcblas
INCLUDE=-I/data/hqchen/software/clibs/include


CODEDIR=/data/hqchen/software/RT1D_PZ/src/

srcs=propagate.c qss.c calc_rates.c integrate.c cosmology.c atomic_data.c spec.c base.c
srcpaths=$(patsubst %.c, $(CODEDIR)/%.c, $(srcs))

srcs2=propagate_flickering.c qss.c calc_rates.c integrate.c cosmology.c atomic_data.c spec.c base.c
srcpaths2=$(patsubst %.c, $(CODEDIR)/%.c, $(srcs2))

run_real_los: 
	$(CC) $(CFLAGS) -o run_real_los run_real_los.c $(srcpaths) $(LIB) $(INCLUDE) -I $(CODEDIR) 
run_with_LC:
	$(CC) $(CFLAGS) -o run_with_LC flickering.c $(srcpaths2) $(LIB) $(INCLUDE) -I $(CODEDIR) 
