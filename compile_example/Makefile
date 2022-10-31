# Compiler
CC=gcc
# Compiler options, including optimization
CFLAGS= -g -Wall -std=c99
LIB= -lm -lgsl
# If your compiler doesn't support -Ofast, replace it with -O3 or -O2.
# (-Ofast is *significantly* faster, so maybe update your compiler!
#  As far as I can tell, it does not substantially impact any of the output.)


CODEDIR=/Users/hqchen/software/RT1D_PZ/src/

srcs=propagate.c qss.c calc_rates.c integrate.c cosmology.c atomic_data.c spec.c base.c
srcpaths=$(patsubst %.c, $(CODEDIR)/%.c, $(srcs))


run_real_los: 
	$(CC) $(CFLAGS) -o run_real_los run_real_los.c $(srcpaths) $(LIB) -I $(CODEDIR) 
