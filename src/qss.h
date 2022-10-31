#ifndef __QSS_H__
#define __QSS_H__

#ifndef MIN
#define MIN(x,y)        (((x) < (y)) ? (x): (y))
#endif
#ifndef MAX
#define MAX(x,y)        (((x) > (y)) ? (x): (y))
#endif
#define cart_alloc(type,size)            (type *)cart_alloc_worker((size)*sizeof(type),__FILE__,__LINE__) 
#define cart_free(ptr)                   { cart_free_worker(ptr,__FILE__,__LINE__); ptr = NULL; }

#include <stdlib.h>

void* cart_alloc_worker(size_t size, const char *file, int line);


typedef struct outArrayBlock_t{
    size_t size;
    size_t current;
    double *t;
    double **y;
}outArrayBlock;



/* Internal storage for intermediate variables */
typedef struct {
	int num_eqn; 
	double *y0;
	double *y1;
	double *rs; 
	double *a0; 
	double *a1;
	double *w0;
	double *w1;
	double *buf;
	void (* rates)(double t, double *y, void *params, double *w, double *a);
	void (* adjust)(double t, double *y, void *params);
} qss_system;

/*
//  Allocated internal storage - call before calling qss_solve
//  Input:
//    num_eqn - number of equations to solve
//    rates - function that returns rates for a given vector of variables y;
//            Equations to be solved are \dot{y} = w - a*y, where w and a may depend on y and t
//            w can be negative, but a must be non-negative; if a is 0, then QSS reduces to second order Runge-Kutta
//            The speedup due to semi-implicit form of the solver is obtained when w and a depend on y slowly.
//            Signature of rates: void (* rates)(double t, double *y, void *params, double *w, double *a);
//            Parameter of rates:
//            t - in: initial time
//            y - in: initial state of variables
//            params - in: optional additional parameters to send, can be NULL
//            w - out: the vector of w values (as in w-a*y)
//            a - out: the vector of a values (incomplete jacobian)
//    adjust - function to correct values after the time step (if needed, for example to make sure some 
//            combination of variables y is conserved exactly), can be NULL
*/
qss_system *qss_alloc( size_t num_eqn,  
		void (* rates)(double, double *, void *, double *, double *),
		void (* adjust)(double, double *, void *) ); 

/*
//  Frees internal storage after the call to qss_solve
*/
void qss_free( qss_system *sys );

/*
//  Main solver.
//  Input:
//    sys - in: previously allocated internal storage
//    t_begin - in: initial time
//    delta_t - in: time interval to integrate over
//    y - in/out: initial vector of variables on input, the solution at t=t_begin+delta_t on the output
//    err - in: vector of precision values for each equation to solve; different equations may have
//          different error requirements if, for example, one of the variables is close to 0
//    params: optional parameters to send to rates(...)
*/
void qss_solve( qss_system *sys, double t_begin, double t_end, double y[], const double err[], void *params );

void qss_solve_save( qss_system *sys, double t_begin, double t_end, double y[], 
		const double err[], void *params, outArrayBlock *outArr);

outArrayBlock *output_alloc( size_t num_variable, size_t initSize);
void free_block(outArrayBlock *oldblock);
void output_alloc_extend (outArrayBlock *oldblock, size_t num_variable);
#endif
