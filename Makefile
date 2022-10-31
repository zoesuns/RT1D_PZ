# Compiler
CC=gcc
# Compiler options, including optimization
CFLAGS= -g -Wall -lgsl -lm -std=c99
# If your compiler doesn't support -Ofast, replace it with -O3 or -O2.
# (-Ofast is *significantly* faster, so maybe update your compiler!
#  As far as I can tell, it does not substantially impact any of the output.)

all: 
	$(CC) $(CFLAGS) -o rtexe propagate.c test_func.c integrate.c qss.c atomic_data.c base.c spec.c calc_rates.c -lgsl -lm

testP: 
	$(CC) $(CFLAGS) -o testP test.c base.c atomic_data.c spec.c calc_rates.c integrate.c  propagate.c qss.c test_func.c -lgsl -lm

qssmake:
	$(CC) $(CFLAGS) -o testqss test_qss.c integrate.c qss.c atomic_data.c base.c spec.c calc_rates.c  test_func.c -lgsl -lm

sorttest:
	$(CC) $(CFLAGS) -o sorttest sort.c qss.c atomic_data.c base.c spec.c calc_rates.c  -lgsl -lm

clean:
	rm rtexe 

