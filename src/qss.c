#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qss.h"
#include "base.h"

#include "base.h"
#include "atomic_data.h"
#include "spec.h"
#include "calc_rates.h"
#include "qss.h"
#include "integrate.h"
#include "cosmology.h"

void* cart_alloc_worker(size_t size, const char *file, int line)
{
  void *ptr;

  if(size > 0)
    {
      ptr = malloc( size );

      if(ptr == NULL)
	{
	  printf( "Failure allocating\n" );
	}

      return ptr;
    }
  else
    {
      return NULL;
    }
}


void* cart_realloc_worker(void *oldptr, size_t size, const char *file, int line)
{
  void *ptr;

  if(size > 0)
    {
      ptr = realloc(oldptr,size);

      if(ptr == NULL)
	{
	  printf( "Failure reallocating\n" );
	}

      return ptr;
    }
  else
    {
      return NULL;
    }
}


void cart_free_worker(void *ptr, const char *file, int line)
{
  if(ptr != NULL)
    {
      free(ptr);
    }
}




outArrayBlock *output_alloc( size_t num_variable, size_t initSize){

    outArrayBlock *block = cart_alloc(outArrayBlock,1);
    block->size=initSize;
    block->current=0;
    block->t = cart_alloc(double, initSize);
    block->y = (double **)cart_alloc_worker((initSize)*sizeof(double*),__FILE__,__LINE__); 
    for (int i=0; i<initSize; i++){
        block->y[i]=cart_alloc(double, num_variable);
    }
    return block;
}

void free_block(outArrayBlock *oldblock){
    for (int i=0; i<oldblock->size; i++){
        cart_free(oldblock->y[i]);
    }
    cart_free(oldblock->t);
    cart_free(oldblock);
}



void output_alloc_extend (outArrayBlock *oldblock, size_t num_variable){
    size_t newSize=2*oldblock->size;

    double *newt=cart_alloc(double, newSize);
    //memcpy(newt, oldblock->t, sizeof(oldblock->t));
    memcpy(newt, oldblock->t, sizeof(double)*oldblock->size);
    double **newY=(double **)cart_alloc_worker((newSize)*sizeof(double*),__FILE__,__LINE__);  
    int i=0;
    for (i=0; i<oldblock->size; i++){
        newY[i]=cart_alloc(double, num_variable);
        memcpy(newY[i], oldblock->y[i], sizeof(double)*num_variable);
        cart_free(oldblock->y[i]);
    }
    for (i=oldblock->size; i<newSize; i++){
        newY[i]=cart_alloc(double, num_variable);

    }


    double *tmpt=oldblock->t;

    oldblock->size=newSize;
    oldblock->t=newt;
    oldblock->y=newY;

    cart_free(tmpt);
}


qss_system *qss_alloc( size_t num_eqn,  
		void (* rates)(double, double *, void *, double *, double *),
		void (* adjust)(double, double *, void *) ) {

	qss_system *sys = cart_alloc(qss_system, 1 );
	sys->buf = cart_alloc(double, 7*num_eqn );
	sys->y0 = sys->buf;
	sys->y1 = &sys->buf[num_eqn];
	sys->rs = &sys->buf[2*num_eqn];
	sys->a0 = &sys->buf[3*num_eqn];
	sys->w0 = &sys->buf[4*num_eqn];
	sys->a1 = &sys->buf[5*num_eqn];
	sys->w1 = &sys->buf[6*num_eqn];
	sys->num_eqn = num_eqn;
	sys->rates = rates;
	sys->adjust = adjust;
	return sys;
}


void qss_free( qss_system *sys ) {
	cart_free( sys->buf );
	cart_free( sys );
}

double qss_alpha( double tau ) {
	if ( tau < 0.01 ) {
		return 0.5 + tau/12 - tau*tau*tau/720;
	} else {
		return (1-(1-exp(-tau))/tau)/(1-exp(-tau));
	}
}

void qs1_step( qss_system *sys, double t, double dt, double yf[], void *params ) {
	int i;
	double tau, alp, pbar;

	for ( i = 0; i < sys->num_eqn; i++ ) {
		tau = dt*sys->a0[i];
		alp = qss_alpha(tau);
		sys->y1[i] = sys->y0[i] +
			dt*(sys->w0[i]-sys->a0[i]*sys->y0[i])/(1+alp*tau);
	}

	/* compute rates for corrector step */
	sys->rates( t, sys->y1, params, sys->w1, sys->a1 );

	for ( i = 0; i < sys->num_eqn; i++ ) {
		pbar = 0.5*(sys->a0[i]+sys->a1[i]);
		tau = pbar*dt;
		alp = qss_alpha(tau);
		yf[i] = sys->y0[i] +
			dt*(alp*sys->w1[i]+(1-alp)*sys->w0[i]- pbar*sys->y0[i]) /
			(1+alp*tau);
		sys->rs[i] = yf[i] - sys->y1[i];
	}
}

void qsn_step( qss_system *sys, double t, double dt, double yf[], void *params ) {
	int i;
	int nIter, nCor;
	double rCor;
	double tau, alp, pbar;

	for ( i = 0; i < sys->num_eqn; i++ ) {
		tau = dt*sys->a0[i];
		alp = qss_alpha(tau);
		sys->y1[i] = sys->y0[i] +
			dt*(sys->w0[i]-sys->a0[i]*sys->y0[i])/(1+alp*tau);
		sys->rs[i] = sys->y1[i];
	}

	nIter = 0;
	nCor = 2;
	rCor = 1.0;

	do {
		/* compute rates for corrector step */
		sys->rates( t, sys->y1, params, sys->w1, sys->a1 );

		for ( i = 0; i < sys->num_eqn; i++ ) {
			pbar = 0.5*(sys->a0[i]+sys->a1[i]);
			tau = pbar*dt;
			alp = qss_alpha(tau);
			yf[i] = sys->y0[i] +
				dt*(alp*sys->w1[i]+(1-alp)*sys->w0[i]-pbar*sys->y0[i]) /
				(1+alp*tau);
		}

		/* Nick's black magic, designed to advance further when mode is oscillating */
		if ( (yf[0]-sys->y0[0])*(sys->rs[0]-sys->y0[0]) < 0.0 ) {
			nCor = 4;
			rCor = 0.2;
		}

		nIter++;

		if ( nIter < nCor ) {
			for ( i = 0; i < sys->num_eqn; i++ ) {
				sys->y1[i] = (1-rCor)*sys->y1[i] + rCor*yf[i];
			}
		}
	} while ( nIter < nCor );

	for ( i = 0; i < sys->num_eqn; i++ ) {
		sys->rs[i] = yf[i] - sys->rs[i];
	}
}

void qss_solve( qss_system *sys, double t_begin, double t_end, double y[], 
		const double err[], void *params ) {
	int i;
	int nstep = 0;
	double errmax;
	double erri;

	double t = t_begin;
	double dt = t_end-t_begin;

	do {
		for ( i = 0; i < sys->num_eqn; i++ ) {
			sys->y0[i] = y[i];
		}

		/* compute initial rates */
		sys->rates( t, sys->y0, params, sys->w0, sys->a0 );

		while ( 1 ) {
			/* take step */
			if ( nstep < 300 ) {
				qs1_step( sys, t, dt, y, params );
			} else {
				qsn_step( sys, t, dt, y, params );
			}

			/* check for errors */
			errmax = fabs(sys->rs[0])/MAX(fabs(y[0]),1e-30)/err[0];
			for ( i = 1; i < sys->num_eqn; i++ ) {
				erri = fabs(sys->rs[i])/MAX(fabs(y[i]),1e-30)/err[i];
				if ( erri > errmax ) {
					errmax = erri;
				}
			}
            //printf("error=%e\n",errmax);
			if ( errmax > 1.0 ) {
				dt = MIN( 0.9*dt/sqrt(2.+errmax), t_end-t );
			} else {
				t += dt;
				nstep++;
				dt = MIN( 0.9*dt/pow(0.01+errmax,0.3), t_end-t );
				break;
			}
		}
		/* adjust variables for max/min */
		if(sys->adjust != NULL) sys->adjust( t, y, params );
	} while ( t < t_end && nstep < 100000 ); 
}


void qss_solve_save( qss_system *sys, double t_begin, double t_end, double y[], 
		const double err[], void *params, outArrayBlock *tmpOut) {
	int i;
	int nstep = 0;
	double errmax;
	double erri;

	double t = t_begin;
	double dt = t_end-t_begin;


	do {
		for ( i = 0; i < sys->num_eqn; i++ ) {
			sys->y0[i] = y[i];
		}

		/* compute initial rates */
		sys->rates( t, sys->y0, params, sys->w0, sys->a0 );

		while ( 1 ) {
			/* take step */
			if ( nstep < 300 ) {
				qs1_step( sys, t, dt, y, params );
			} else {
				qsn_step( sys, t, dt, y, params );
			}

			/* check for errors */
			errmax = fabs(sys->rs[0])/MAX(fabs(y[0]),1e-30)/err[0];
			for ( i = 1; i < sys->num_eqn; i++ ) {
				erri = fabs(sys->rs[i])/MAX(fabs(y[i]),1e-30)/err[i];
				if ( erri > errmax ) {
					errmax = erri;
				}
			}




			if ( errmax > 1.0 ) {
				dt = MIN( 0.9*dt/sqrt(2.+errmax), t_end-t );
			} else {
				t += dt;
                ////////////// save //////////
                if (tmpOut->current>=tmpOut->size){
                    output_alloc_extend (tmpOut, sys->num_eqn);
                }


                tmpOut->t[tmpOut->current]=t;
                for (i=0; i<sys->num_eqn; i++){
                    tmpOut->y[tmpOut->current][i]=y[i];
                }
                tmpOut->current+=1;
                //////////////////////////////
				nstep++;
				dt = MIN( 0.9*dt/pow(0.01+errmax,0.3), t_end-t );
				break;
			}
		}
		/* adjust variables for max/min */
		if(sys->adjust != NULL) sys->adjust( t, y, params );
	} while ( t < t_end && nstep < 100000 ); 
}

