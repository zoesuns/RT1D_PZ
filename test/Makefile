# Compiler
CC=gcc
# Compiler options, including optimization
CFLAGS= -g -Wall -std=c99
LIB= -lm -lgsl
# If your compiler doesn't support -Ofast, replace it with -O3 or -O2.
# (-Ofast is *significantly* faster, so maybe update your compiler!
#  As far as I can tell, it does not substantially impact any of the output.)

CODEDIR="../"

srcs=propagate.c qss.c calc_rates.c integrate.c cosmology.c atomic_data.c spec.c base.c
srcpaths=$(patsubst %.c, $(CODEDIR)/%.c, $(srcs))


IO: 
	$(CC) $(CFLAGS) -o test_IO test_IO.c ../atomic_data.c ../base.c $(LIB)

atomic_data: 
	$(CC) $(CFLAGS) -o print_atomic_data print_atomic_data.c fake_cells.c ../cosmology.c  ../atomic_data.c ../base.c $(LIB)

cosmology: 
	$(CC) $(CFLAGS) -o print_cosmology print_cosmology.c ../cosmology.c  ../base.c $(LIB)

spectra: 
	$(CC) $(CFLAGS) -o print_spectra print_spectra.c ../atomic_data.c ../spec.c ../base.c $(LIB)

photoGH: 
	$(CC) $(CFLAGS) -o test_photoGH test_photoGH.c fake_cells.c ../qss.c ../calc_rates.c ../integrate.c ../atomic_data.c ../spec.c ../base.c $(LIB)

integrate: 
	$(CC) $(CFLAGS) -o test_integrate test_integrate.c fake_cells.c ../qss.c ../calc_rates.c ../integrate.c ../atomic_data.c ../spec.c ../base.c $(LIB)

strange_cell: 
	$(CC) $(CFLAGS) -o strange_cell strange_cell.c fake_cells.c ../qss.c ../calc_rates.c ../integrate.c ../cosmology.c ../atomic_data.c ../spec.c ../base.c $(LIB)

propagate: 
	$(CC) $(CFLAGS) -o test_propagate test_propagate.c fake_cells.c ../propagate.c ../qss.c ../calc_rates.c ../integrate.c ../cosmology.c ../atomic_data.c ../spec.c ../base.c $(LIB)

fast_propagate: 
	$(CC) $(CFLAGS) -o test_fast_propagate test_fast_propagate.c fake_cells.c ../fast_propagate.c ../qss.c ../calc_rates.c ../integrate.c ../cosmology.c ../atomic_data.c ../spec.c ../base.c $(LIB)

los: 
	$(CC) $(CFLAGS) -o test_los test_los.c fake_cells.c ../propagate.c ../qss.c ../calc_rates.c ../integrate.c ../atomic_data.c ../cosmology.c ../spec.c ../base.c $(LIB)

HII_sphere: 
	$(CC) $(CFLAGS) -o pureHII test_HII_sphere.c fake_cells.c $(srcpaths) $(LIB) -I $(CODEDIR) 

CGM1: 
	$(CC) $(CFLAGS) -o CGM_02_04pMpc test_CGM_Ifront1.c fake_cells.c $(srcpaths) $(LIB) -I $(CODEDIR) 

CGM2: 
	$(CC) $(CFLAGS) -o CGM_04_05pMpc test_CGM_Ifront2.c fake_cells.c $(srcpaths) $(LIB) -I $(CODEDIR) 




