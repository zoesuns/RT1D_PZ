#include "config.h"
#ifdef HYDRO 

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cooling.h"
#include "cosmology.h"
#include "hydro.h"
#include "iterators.h"
#include "qss.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "hydro_step.h"
#include "hydro_sgst_step.h"
#include "rt_step.h"
#include "step.h"
#include "plugin.h"


#ifndef HYDRO_CHUNK_SIZE
#define HYDRO_CHUNK_SIZE        65536
#endif /* HYDRO_CHUNK_SIZE */

extern int pressure_floor_min_level; /* NG: that used to be MinL_Jeans define */
extern float pressure_floor_factor;
float pressure_floor;

extern int pressureless_fluid_eos;           /* NG: that used to be PRESSURELESS_FLUID define */
extern int apply_lapidus_viscosity;          /* NG: that used to be LAPIDUS define */
extern int smooth_density_gradients;         /* NG: that used to be DENSGRADSMOOTH define */

extern float gas_density_floor;
extern float gas_temperature_floor;          /* NG: that used to be T_min define */

#ifdef BLASTWAVE_FEEDBACK
extern double blastwave_time_floor; 
extern double blastwave_time_cut;
#endif /* BLASTWAVE_FEEDBACK */

extern float fixed_metallicity_log;
extern double shielding_min_density_log;
extern double shielding_jeans_temperature_ceiling;

#if defined(COSMOLOGY) && defined(REFINEMENT)
extern double fixed_proper_resolution;
#endif /* COSMOLOGY && REFINEMENT */


#define backup_hvar(c,v)	(backup_hvars[c][v])

float backup_hvars[num_cells][num_hydro_vars-2];
float ref[num_cells];
int backup_dirty[num_cells];
#ifdef HYDRO_TRACERS_MONTE_CARLO
float remaining_density[num_cells];
#else /* HYDRO_TRACERS_MONTE_CARLO */
float remaining_density; // dummy variable for #pragma omp parallel
#endif /* HYDRO_TRACERS_MONTE_CARLO */

#ifdef GRAVITY_IN_RIEMANN
void fluxh( double dtx, double dtx2, double v[num_state_vars][4], double g[2], double c[2], double f[num_hydro_vars-1] );
#else
void fluxh( double dtx, double dtx2, double v[num_state_vars][4], double c[2], double f[num_hydro_vars-1] );
#endif

void lapidus( double dtx2, int L1, int R1, int sweep_direction, int mj3, int mj4, int mj5, double v[num_state_vars][4], double f[num_hydro_vars-1] );

const int sweep_dir[2][nDim] = { { 1, 3, 5 }, { 5, 3, 1 } };

int sweep_direction;
int sweep_dimension;

int j3,j4,j5;
int mj3, mj4, mj5;
double dtx;
double dtx2;
double dxi;
double dxi2;

int momentum_permute[2*nDim][nDim] = {  
	{ 0, 1, 2 }, { 0, 1, 2 },
	{ 1, 0, 2 }, { 1, 0, 2 },
	{ 2, 1, 0 }, { 2, 1, 0 } };


void apply_hydro_fluxes( int icell, double factor, double dxi_factor, double f[ /* num_hydro_vars-1 */ ] );
void hydro_sweep_1d( int level );
#ifdef GRAVITY
void hydro_apply_gravity( int level );
#endif /* GRAVITY */
void compute_hydro_fluxes( int cell_list[4], double f[ /* num_hydro_vars-1 */ ] );
void hydro_advance_internalenergy(int level);

void hydro_step( int level ) {
	int dir;

	if ( pressure_floor_min_level >= 0 ) {
		/* artificial pressure floor */
		double len = 0.0;
 
#if defined(COSMOLOGY) && defined(REFINEMENT)
                if ( fixed_proper_resolution > 0.0 ) {
			len = fixed_proper_resolution*constants->kpc/units->length;
			if( len < cell_size[level] ) len = (level >= pressure_floor_min_level) ? cell_size[level] : 0.0;
                } else {
#endif
			if ( level >= pressure_floor_min_level ) len = cell_size[level];
#if defined(COSMOLOGY) && defined(REFINEMENT)
                }
#endif
		pressure_floor = pressure_floor_factor * constants->G*pow(units->density*units->length,2.0)/units->energy_density * len *len;

		start_time( PLUGIN_TIMER );
		PLUGIN_POINT(PressureFloor,(level, len, &pressure_floor));
		end_time( PLUGIN_TIMER );
	} else {
		pressure_floor = 0.0;
	}

	dtx = dtl[level] * cell_size_inverse[level];
	dxi = cell_size_inverse[level];
	dxi2 = 0.5*cell_size_inverse[level];
	dtx2 = 0.5*dtx;

	start_time( HYDRO_TIMER );

	for ( dir = 0; dir < nDim; dir++ ) {
		start_time( WORK_TIMER );

		sweep_direction = sweep_dir[level_sweep_dir[level]][dir];
		sweep_dimension = (sweep_direction-1)/2;

		j3 = momentum_permute[sweep_direction][0];
		j4 = momentum_permute[sweep_direction][1];
		j5 = momentum_permute[sweep_direction][2];

		mj3 = j3+2;
		mj4 = j4+2;
		mj5 = j5+2;

		/* compute fluxes across cell interfaces long sweep_dimension */
		hydro_sweep_1d( level );

#if defined(GRAVITY) && (!defined(GRAVITY_IN_RIEMANN))
		hydro_apply_gravity( level );
#endif /* GRAVITY && !GRAVITY_IN_RIEMANN */

		end_time( WORK_TIMER );

		/*
		//  Dry-run past the CFL violation to avoid the solution blowing up
		*/
		if (!cfl_violation_occurred) {
			if ( dir == nDim - 1 ) {
				hydro_copy_vars( level, HYDRO_RESTORE_ALL );

#ifdef SGS_TURBULENCE
				hydro_sgst_advance_turbulent_energy( level );
#endif /* SGS_TURBULENCE */

				hydro_advance_internalenergy( level );
				hydro_magic( level );
				hydro_eos( level );

#ifdef SGST_SHEAR_IMPROVED
				hydro_sgst_update_mean_velocities( level );
#endif /* SGST_SHEAR_IMPROVED */

				hydro_split_update( level );
			} else {
				hydro_copy_vars( level, HYDRO_RESTORE_CLEAN );
			    hydro_magic( level );
			    hydro_eos( level );
			}
		}

		start_time( HYDRO_UPDATE_TIMER );
		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
		end_time( HYDRO_UPDATE_TIMER );
	}

	/* update sweep direction */
	level_sweep_dir[level] = (level_sweep_dir[level]+1)%2;

	end_time( HYDRO_TIMER );
}

void hydro_sweep_1d( int level ) {
	int i, j;
	int icell;
	int L2, L1, R1, R2;
	int count;
	int num_level_cells;
	int *level_cells;

	int cell_list[HYDRO_CHUNK_SIZE][4];
	double f[HYDRO_CHUNK_SIZE][num_hydro_vars-1];

	count = 0;

	if ( level == min_level ) {
		select_level( min_level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
	} else {
		select_level( level, CELL_TYPE_ANY_LEAF, &num_level_cells, &level_cells );
	}

	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		/* calculate neighbors */
		L1 = cell_neighbor( icell, reverse_direction[sweep_direction]);
		R1 = cell_neighbor( icell, sweep_direction );

		if ( ( cell_is_local(icell) && cell_is_leaf(R1) ) ||
				(!cell_is_local(icell) && R1 != -1 && cell_is_local(R1) && cell_is_leaf(R1) ) ) {

			R2 = cell_neighbor( R1, sweep_direction );

			cell_list[count][0] = L1;
			cell_list[count][1] = icell;
			cell_list[count][2] = R1;
			cell_list[count][3] = R2;
			count++;
		}

		if ( ( level == min_level && !cell_is_local(L1) && cell_is_leaf(L1) ) || 
				( L1 != -1 && cell_level(L1) == level - 1 && 
				( cell_is_local(icell) || cell_is_local(L1) ) ) ) {

			L2 = cell_neighbor( L1, reverse_direction[sweep_direction] );

			cell_list[count][0] = L2;
			cell_list[count][1] = L1;
			cell_list[count][2] = icell;
			cell_list[count][3] = R1;
			count++;
		}

		/* must be HYDRO_CHUNK_SIZE-1 since each cell can add 2 interfaces */
		if ( i == num_level_cells-1 || count >= HYDRO_CHUNK_SIZE-1 ) {
#pragma omp parallel for default(none), private(j), shared(cell_list,f,count,dtx,dtx2,sweep_direction,j3,j4,j5)
			for ( j = 0; j < count; j++ ) {
				compute_hydro_fluxes( cell_list[j], f[j] );
			}

 			/* apply fluxes to left cells */
#pragma omp parallel for default(none), private(j,icell), shared(cell_list,count,level,f,dxi,dxi2)
			for ( j = 0; j < count; j++ ) {
				icell = cell_list[j][1];

				if ( cell_is_local(icell) && cell_level(icell) == level ) {
					apply_hydro_fluxes( icell, -1.0, dxi, f[j] );
				}
			}

#pragma omp parallel for default(none), private(j,icell), shared(count,cell_list,level,f,dxi,dxi2)
			for ( j = 0; j < count; j++ ) {
				icell = cell_list[j][2];

				if ( cell_is_local(icell) && cell_level(icell) == level ) {
					apply_hydro_fluxes( icell, 1.0, dxi, f[j] );
				}
			}

			/* Apply fluxes to higher level cells on both left and right interfaces 
			 * MUST execute serially since lower level cells will appear multiple times! */
			for ( j = 0; j < count; j++ ) {
				icell = cell_list[j][1];

				if ( cell_is_local(icell) && cell_level(icell) < level ) {
					apply_hydro_fluxes( icell, -0.125, dxi2, f[j] );
					backup_dirty[icell] = 1;
				}

				icell = cell_list[j][2];
				if ( cell_is_local(icell) && cell_level(icell) < level ) {
					apply_hydro_fluxes( icell, 0.125, dxi2, f[j] );
					backup_dirty[icell] = 1;
				}

#ifdef HYDRO_TRACERS_MONTE_CARLO
				transfer_hydro_tracers( cell_list[j][1], cell_list[j][2], f[j][0] );
#endif /* HYDRO_TRACERS_MONTE_CARLO */
			}

			count = 0;
		}
	}

	cart_free( level_cells );
}

#if defined(GRAVITY) && (!defined(GRAVITY_IN_RIEMANN))
void hydro_apply_gravity( int level ) {
	int i, icell;
	int num_level_cells;
	int *level_cells;
	double gravadd;

	/* now we need to apply a gravity correction */
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(icell,gravadd), shared(num_level_cells,level_cells,cell_child_oct,backup_hvars,cell_vars_data,sweep_dimension,mj3)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		gravadd = backup_hvar(icell,0) * cell_accel(icell,sweep_dimension);
		backup_hvar(icell,1) += cell_accel(icell,sweep_dimension) *
			( backup_hvar(icell,mj3) + 0.5 * gravadd );
		backup_hvar(icell,mj3) += gravadd;
	}

	cart_free( level_cells );
}
#endif /* GRAVITY && !GRAVITY_IN_RIEMANN */

void hydro_eos( int level ) {
    int i,j;
	int icell;
	int num_level_cells;
	int *level_cells;
	double kinetic_energy;
	/*
	//  Dereference for efficiency
	*/
#ifdef STAR_FORMATION
	float (*extra_pressure)(int cell) = sf_feedback_particle->extra_pressure;
#else
	float (*extra_pressure)(int cell) = NULL;
#endif

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,j,icell,kinetic_energy), shared(num_level_cells,level_cells,cell_child_oct,cell_vars_data,constants,pressureless_fluid_eos,extra_pressure)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		kinetic_energy = cell_gas_kinetic_energy(icell);

		if(pressureless_fluid_eos)
		  {
		    cell_gas_pressure(icell) = 1e-20;
		    cell_gas_internal_energy(icell) = cell_gas_pressure(icell) / (constants->gamma-1.0);
		    cell_gas_energy(icell) = cell_gas_internal_energy(icell) + kinetic_energy;
		  }
		else
		  {
		    cell_gas_internal_energy(icell) = MAX( cell_gas_internal_energy(icell), 0.0 );
		    for ( j = 0; j < num_extra_energy_variables; j++ ) {
		    	cell_extra_energy_variables(icell,j) = MAX( cell_extra_energy_variables(icell,j), 0.0 );
		    }
		    cell_gas_energy(icell) = MAX( kinetic_energy, cell_gas_energy(icell) );
		    cell_gas_pressure(icell) = MAX( (cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell), 0.0 );

		    if(extra_pressure != NULL) cell_gas_pressure(icell) += extra_pressure(icell);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		    cell_electron_internal_energy(icell) = MIN( cell_electron_internal_energy(icell), cell_gas_internal_energy(icell)*constants->wmu/constants->wmu_e );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
		  }
	}

	cart_free( level_cells );

	end_time( WORK_TIMER );
}


#if defined(COOLING) && !defined(RADIATIVE_TRANSFER)

#ifndef OLDSTYLE_COOLING_EXPLICIT_SOLVER
void qss_getcooling ( double t, double *y, void *params, double *w, double *a) {
	double nHlog = ((double *)params)[0];
	double Tfac_cell = ((double *)params)[1];
	double Zlog = ((double *)params)[2];
	double rhog2 = ((double *)params)[3];
	double unit_cl = ((double *)params)[6];
	double t0 = ((double *)params)[7];
	double Hdum = ((double *)params)[8];
	double f_curr = 1 + Hdum*(t-t0);
	double etmp = rhog2*f_curr*unit_cl;
	double nTcode = y[0]/(f_curr*f_curr);

	cooling_t coolrate = cooling_rate(nHlog,Tfac_cell*nTcode,Zlog);
	
	a[0] = etmp*coolrate.Cooling/ y[0];
	w[0] = etmp*coolrate.Heating;

#ifdef COOLING_DUST_SHIELDING
	w[0] *= exp( MAX( ((double *)params)[11]*sqrt(nTcode), ((double *)params)[12] ) );
#endif /* COOLING_DUST_SHIELDING */

	a[0] *= ((double *)params)[9];
	w[0] *= ((double *)params)[10];
}

void adjust_internalenergy( double t, double *y, void *params ) {
  /* RL: put temperature/internal energy floor in here??? */
	double Emin_cell = ((double *)params)[5];
	if (y[0] < Emin_cell) y[0] = Emin_cell;
}

void hydro_apply_cooling(int level, int num_level_cells, int *level_cells) {
	int i, j;
	int icell;
	double t_begin, t_stop;
	double Z, Hdum;
	double Tfac, Tfac_cell, Emin_cell;
	double log_fdust_fac, log_fdust_max_fac, nHfac;
#ifdef BLASTWAVE_FEEDBACK
	double blastwave_time;
#endif /* BLASTWAVE_FEEDBACK */
	double Eminfac;
	double rhog2, nHlog, nH;
	double e_curr;
	double unit_cl = units->time*pow(constants->XH*units->number_density,2.0)/units->energy_density;
	double cooling_multiplier, heating_multiplier;
	double err[1] = { 1e-2 };
	double params[13];
	qss_system *sys;

	t_begin	= tl[level];
	t_stop = tl[level] + dtl[level];

#ifdef COSMOLOGY
	Hdum = (abox_from_tcode(t_stop)/abox[level] - 1) / dtl[level];
#else
	Hdum = 0.0;
#endif

	/* Note: removed 10^4 K normalization since it wasn't used - DHR */
	Tfac = units->temperature*constants->wmu*( constants->gamma-1 )/constants->K;
	Eminfac = gas_temperature_floor/(units->temperature*constants->wmu*(constants->gamma-1));
	nHfac = constants->XH*units->number_density/constants->cc;

	/* eqs. 15, 16 & 27 in Safranek-Shrader et al. (2016) with gamma = 2.5 and sigma_{d,V} = 5.3e-22 cm^2: */
	log_fdust_fac = -2.5*5.3e-22*nHfac*units->velocity*sqrt(M_PI*constants->gamma*(constants->gamma-1)/(constants->G*units->density));
	log_fdust_max_fac = log_fdust_fac*sqrt((shielding_jeans_temperature_ceiling*constants->K)/(units->temperature*constants->wmu*(constants->gamma-1)));


#ifdef BLASTWAVE_FEEDBACK
#pragma omp parallel default(none), shared(num_level_cells,level_cells,level,t_begin,Tfac,Eminfac,nHfac,log_fdust_fac,log_fdust_max_fac,shielding_min_density_log,units,t_stop,cell_child_oct,err,constants,cell_vars_data,Hdum,unit_cl,fixed_metallicity,blastwave_time_cut,blastwave_time_floor,plugins), private(i,j,icell,rhog2,nHlog,Z,Tfac_cell,e_curr,Emin_cell,params,sys,blastwave_time,cooling_multiplier,heating_multiplier)
#else
#pragma omp parallel default(none), shared(num_level_cells,level_cells,level,t_begin,Tfac,Eminfac,nHfac,log_fdust_fac,log_fdust_max_fac,shielding_min_density_log,units,t_stop,cell_child_oct,err,constants,cell_vars_data,Hdum,unit_cl,fixed_metallicity,plugins), private(i,j,icell,rhog2,nHlog,Z,Tfac_cell,e_curr,Emin_cell,params,sys,cooling_multiplier,heating_multiplier)
#endif /* BLASTWAVE_FEEDBACK*/
	{
	  sys = qss_alloc( 1, &qss_getcooling, &adjust_internalenergy );

#pragma omp for
	  for ( i = 0; i < num_level_cells; i++ ) {
		  icell = level_cells[i];
		  if ( cell_is_leaf(icell) ) {
#ifdef BLASTWAVE_FEEDBACK
		    blastwave_time = cell_blastwave_time(icell) / cell_gas_density(icell);
		    if(blastwave_time <= blastwave_time_cut){ 
#endif /* BLASTWAVE_FEEDBACK */

		    cell_gas_gamma(icell) = constants->gamma;
		    rhog2 = cell_gas_density(icell)*cell_gas_density(icell);
		    /* take code density -> log10(n_H [cm^-3]) */

		    nHlog = log10(nHfac*cell_gas_density(icell));

#ifdef ENRICHMENT
		    Z = MAX(1.0e-10,cell_gas_metal_density(icell)/(constants->Zsun*cell_gas_density(icell)));
#else
		    Z = fixed_metallicity;
#endif /* ENRICHMENT */
			  
		    Tfac_cell = Tfac/cell_gas_density(icell);
		    Emin_cell = Eminfac*cell_gas_density(icell);
		    
		    e_curr = cell_gas_internal_energy(icell);

		    cooling_multiplier = heating_multiplier = 1.0;
		    PLUGIN_POINT(ModifyCooling,(icell, level, &cooling_multiplier, &heating_multiplier));

		    params[0] = nHlog; //
		    params[1] = Tfac_cell; // to get the cooling rate...
		    params[2] = log10(Z); //
		    params[3] = rhog2;
		    params[4] = e_curr;
		    params[5] = Emin_cell;
		    params[6] = unit_cl;
		    params[7] = t_begin;
		    params[8] = Hdum;
		    params[9] = cooling_multiplier;
		    params[10] = heating_multiplier;
#ifdef COOLING_DUST_SHIELDING
		    params[11] = (nHlog > shielding_min_density_log) ? log_fdust_fac*Z : 0.0;
		    params[12] = log_fdust_max_fac*Z*sqrt(cell_gas_density(icell));
#endif /* COOLING_DUST_SHIELDING */
	
		    qss_solve( sys, t_begin, t_stop, &e_curr, err, &params );

#ifdef ADVECT_EXTRA_ENERGIES
		    e_curr = (cell_gas_internal_energy(icell) - MAX(Emin_cell,e_curr));
		    cell_gas_internal_energy(icell) -= e_curr;
		    cell_gas_energy(icell) -= e_curr;
#else /* ADVECT_EXTRA_ENERGIES */
		    cell_gas_internal_energy(icell) = MAX(Emin_cell,e_curr);
		    cell_gas_energy(icell) = cell_gas_kinetic_energy(icell) + cell_gas_internal_energy(icell);
		    for ( j = 0; j < num_extra_energy_variables; j++ ) {
		    	cell_extra_energy_variables(icell,j) = MAX( cell_extra_energy_variables(icell,j), 0.0 );
		    	cell_gas_energy(icell) += cell_extra_energy_variables(icell,j);
		    }
#endif /* ADVECT_EXTRA_ENERGIES */

#ifdef BLASTWAVE_FEEDBACK
			} else { 
				blastwave_time -= dtl[level]*units->time/constants->yr; 
				if(blastwave_time < blastwave_time_cut ){
					blastwave_time = blastwave_time_floor;
				}
				cell_blastwave_time(icell) = cell_gas_density(icell) * blastwave_time;
			}
#endif /* BLASTWAVE_FEEDBACK */
		  }
		}

		qss_free(sys);
	} /* END omp parallel block */
}

#else /* OLDSTYLE_COOLING_EXPLICIT_SOLVER */

void hydro_cool_one_cell(int icell, double t_begin, double t_stop, double Hdum, double Zlog, double nHlog, double rhog2, double Tfact_cell, double Emin_cell, double unit_cl) {
	int continue_cooling;
	double t_curr, f_curr;
	double dE;
	double dt_e, ei1, T_gas;

#define dstep	(0.01)

	continue_cooling = 1;
	t_curr = t_begin;

	/* integrate cooling using smaller timestep */
	while ( ( t_curr < t_stop ) && continue_cooling ) {
		f_curr = 1 + Hdum*(t_curr-t_begin);
		T_gas = Tfact_cell * cell_gas_internal_energy(icell) / ( f_curr * f_curr );

		/* compute new timestep */
		dE = unit_cl*cooling_rate( nHlog, T_gas, Zlog );
		dE *= -rhog2 * f_curr;
		dt_e = MIN( dstep * fabs( cell_gas_internal_energy(icell) / dE ), t_stop - t_curr );

		ei1 = MAX( cell_gas_internal_energy(icell) + 0.5 * dE * dt_e, Emin_cell );
		T_gas = Tfact_cell * ei1 / ( f_curr * f_curr );

		dE = unit_cl*cooling_rate( nHlog, T_gas, Zlog );
		dE *= -rhog2 * f_curr * dt_e;
		/* adjust cell energies */
		cell_gas_internal_energy(icell) += dE;
		cell_gas_energy(icell) += dE;

		/* stop if we hit energy minimum */
		if ( cell_gas_internal_energy(icell) < Emin_cell ) {
			continue_cooling = 0;
		}

		cell_gas_internal_energy(icell) = MAX( Emin_cell, cell_gas_internal_energy(icell) );
		cell_gas_energy(icell) = MAX( Emin_cell, cell_gas_energy(icell) );

		/* advance timestep */
		t_curr += dt_e;
	}
}

void hydro_apply_cooling(int level, int num_level_cells, int *level_cells) {
	int i;
	int icell;
	double t_begin, t_stop;
	double Zlog, Hdum;
	double Tfac, Tfac_cell, Emin_cell;
	double Eminfac;
	double rhog2, nHlog;
	double unit_cl = units->time*pow(constants->XH*units->number_density,2.0)/units->energy_density;

	t_begin	= tl[level];
	t_stop = tl[level] + dtl[level];

#ifdef COSMOLOGY
	Hdum = (abox_from_tcode(t_stop)/abox[level] - 1) / dtl[level];
#else
	Hdum = 0.0;
#endif

	/* Note: removed 10^4 K term since it isn't used - DHR */
	Tfac = units->temperature*constants->wmu*(constants->gamma-1)/constants->K;
	Eminfac = gas_temperature_floor/(units->temperature*constants->wmu*(constants->gamma-1));

#ifdef BLASTWAVE_FEEDBACK
#pragma omp parallel for default(none), private(icell,i,rhog2,nHlog,Zlog,Tfac_cell,Emin_cell,blastwave_time), shared(num_level_cells,level_cells,t_begin,t_stop,Tfac,Eminfac,units,constants,cell_child_oct,cell_vars_data,Hdum,unit_cl,fixed_metallicity_log,blastwave_time_cut,blastwave_time_floor)
#else
#pragma omp parallel for default(none), private(icell,i,rhog2,nHlog,Zlog,Tfac_cell,Emin_cell), shared(num_level_cells,level_cells,t_begin,t_stop,Tfac,Eminfac,units,constants,cell_child_oct,cell_vars_data,Hdum,unit_cl,fixed_metallicity_log)
#endif /* BLASTWAVE_FEEDBACK*/ 
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		if ( cell_is_leaf(icell) ) {
#ifdef BLASTWAVE_FEEDBACK
		  blastwave_time = cell_blastwave_time(icell) / cell_gas_density(icell);
		  if(blastwave_time <= blastwave_time_cut){
#endif /* BLASTWAVE_FEEDBACK */
		
			cell_gas_gamma(icell) = constants->gamma;
			rhog2 = cell_gas_density(icell)*cell_gas_density(icell);
			/* take code density -> log10(n_H [cm^-3]) */
			nHlog = log10(constants->XH*units->number_density*cell_gas_density(icell)/constants->cc);

#ifdef ENRICHMENT
			Zlog = log10(MAX(1.0e-10,cell_gas_metal_density(icell)/(constants->Zsun*cell_gas_density(icell))));
#else
			Zlog = fixed_metallicity_log;
#endif /* ENRICHMENT */

			Tfac_cell = Tfac/cell_gas_density(icell);
			Emin_cell = Eminfac*cell_gas_density(icell);

			hydro_cool_one_cell(icell,t_begin,t_stop,Hdum,Zlog,nHlog,rhog2,Tfac_cell,Emin_cell,unit_cl);
			
#ifdef BLASTWAVE_FEEDBACK
		  }else { 
		    blastwave_time -= dtl[level]*units->time/constants->yr; 
		    if(blastwave_time < blastwave_time_cut ){
		      blastwave_time = blastwave_time_floor;
		    }
		    cell_blastwave_time(icell) = cell_gas_density(icell) * blastwave_time;
		  }
#endif /* BLASTWAVE_FEEDBACK */
		}
	}
}
#endif /* OLDSTYLE_COOLING_EXPLICIT_SOLVER */

#endif /* COOLING && !RADIATIVE_TRANSFER */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
void heating_rates ( double t, double *y, void *params, double *w, double *a) {
	double dEfact = ((double *)params)[0];
	double e_equil = ((double *)params)[1];
	double e_init = ((double *)params)[2];
	double t0 = ((double *)params)[3];
	double Hdum = ((double *)params)[4];
	double f_curr = 1 + Hdum*(t-t0);
	double e_curr = MAX( e_init, y[0] );

	a[0] = dEfact*f_curr*f_curr*pow(e_curr,-1.5);
	w[0] = a[0]*e_equil;
}

void adjust_temperatures( double t, double *y, void *params ) {
	double e_equil = ((double *)params)[1];
	double e_init = ((double *)params)[2];

	if ( y[0] < e_init ) {
		y[0] = e_init;
	} else if ( y[0] > e_equil ) {
		y[0] = e_equil;
	}
}

void hydro_apply_electron_heating(int level, int num_level_cells, int *level_cells) {
	int i;
	int icell;
	double t_begin, t_stop, Hdum;
	double e_equil, e_curr;
	double nfact, Tefact, dEfact;
	double n_5, Te7, dEfact_cell;
	double logcoulomb;

	double err[1] = { 1e-2 };
	double params[5];
	qss_system *sys;

	t_begin	= tl[level];
	t_stop = tl[level] + dtl[level];

#ifdef COSMOLOGY
	Hdum = (abox_from_tcode(t_stop)/abox[level] - 1) / dtl[level];
#else
	Hdum = 0.0;
#endif

	nfact = 1.0e5*units->number_density*(1.0/constants->wmu - 1.0/constants->wmu_e)/constants->cc;
	Tefact = units->temperature*constants->wmu_e*(constants->gamma-1)/1.0e7/constants->K;
	dEfact = pow(Tefact,-1.5)*units->time/(constants->yr*6.3e8)/40.0/(1.0-constants->wmu/constants->wmu_e)/constants->erg; 

#pragma omp parallel default(none), shared(nfact,Tefact,dEfact,Hdum,cell_vars_data,t_begin,t_stop,num_level_cells,level_cells,cell_child_oct,err,constants), private(i,icell,e_equil,e_curr,Te7,n_5,logcoulomb,dEfact_cell,params,sys)
	{
		sys = qss_alloc( 1, &heating_rates, &adjust_temperatures );

#pragma omp for
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf(icell) ) {
				e_equil = cell_gas_internal_energy(icell)*constants->wmu/constants->wmu_e;
				e_curr = cell_electron_internal_energy(icell);
				Te7 = Tefact * cell_electron_internal_energy(icell) / cell_gas_density(icell); /* a^2 Te/10^7 K */
		
				n_5 = nfact * cell_gas_density(icell); 
				logcoulomb = MAX( 30.0, 37.8 + log(Te7) - 0.5*log(n_5) );
				dEfact_cell = dEfact*n_5*logcoulomb*pow( cell_gas_density(icell), 1.5);
		
				params[0] = dEfact_cell;
				params[1] = e_equil;
				params[2] = e_curr;
				params[3] = t_begin;
				params[4] = Hdum;
	
				qss_solve( sys, t_begin, t_stop, &e_curr, err, &params );
	
				cell_electron_internal_energy(icell) = MIN( e_curr, e_equil );
			}
		}

		qss_free(sys);
	} /* END omp parallel block */
}
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef EXTRA_PRESSURE_SOURCE
void hydro_zero_extra_source_vars(int level, int num_level_cells, int *level_cells) {
    int i,j, icell; 
#pragma omp parallel for default(none), shared(level,num_level_cells,level_cells,cell_child_oct,cell_vars_data,dtl), private(i,j,icell)
    for ( i = 0; i < num_level_cells; i++ ) {
	icell = level_cells[i];
	for ( j = 0; j < num_extra_source_vars; j++ ) {
	    cell_extra_source_variables(icell,j) = 0.0;
	}
    }
}
#endif /* EXTRA_PRESSURE_SOURCE */

void hydro_advance_internalenergy( int level ) {
	int i,j;
	int icell;
	int num_level_cells;
	int *level_cells;
	double kinetic_energy;
	double energy;
	double gamma1, div, div_dt;

#if num_extra_energy_variables > 0
	double internal_energy_factor;
#endif /* num_extra_energy_variables > 0 */

	start_time( WORK_TIMER );

	div_dt = dtl[level] / 3.0;

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#if num_extra_energy_variables > 0
#pragma omp parallel for default(none), private(icell,i,j,kinetic_energy,energy,internal_energy_factor,gamma1,div), shared(num_level_cells,level_cells,cell_child_oct,cell_vars_data,ref,div_dt,extra_energy_gammas)
#else
#pragma omp parallel for default(none), private(icell,i,j,kinetic_energy,energy,gamma1,div), shared(num_level_cells,level_cells,cell_child_oct,cell_vars_data,ref,div_dt,extra_energy_gammas)
#endif /* num_extra_energy_variables > 0 */
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		/* P dV term */
		gamma1 = cell_gas_gamma(icell) - 1.0;
		div = 1.0 + gamma1 * ref[icell] * div_dt;
		cell_gas_internal_energy(icell) = MAX( 1.0e-30, cell_gas_internal_energy(icell)*div*div*div);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		cell_electron_internal_energy(icell) = MAX( 1.0e-30, cell_electron_internal_energy(icell)*div*div*div );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
		for ( j = 0; j < num_extra_energy_variables; j++ ) {
			gamma1 = extra_energy_gamma(j) - 1.0; 
			div = 1.0 + gamma1 * ref[icell] * div_dt;
			cell_extra_energy_variables(icell,j) = MAX( 0.0, cell_extra_energy_variables(icell,j)*div*div*div );
		}

		/* synchronize internal and total energy */
#ifndef ADVECT_EXTRA_ENERGIES
		kinetic_energy = cell_gas_kinetic_energy(icell);
		energy = cell_gas_energy(icell);

		/* we trust energy over internal energy since it's computed using
		 * the riemann solver rather than just advection equation, so if
		 * internal energy is sufficiently large then compute it from 
		 * e = E - rho * v**2 /2 */
		if ( ( energy - kinetic_energy ) / energy > 1e-3 ) {
#if num_extra_energy_variables > 0
			/* Eint and extra energy variables are rescaled proportionally */
			internal_energy_factor = cell_gas_internal_energy(icell);
			for ( j = 0; j < num_extra_energy_variables; j++ ) {
				internal_energy_factor += cell_extra_energy_variables(icell,j);
			}
			internal_energy_factor = ( energy - kinetic_energy ) / internal_energy_factor;

			cell_gas_internal_energy(icell) *= internal_energy_factor;
			for ( j = 0; j < num_extra_energy_variables; j++ ) {
				cell_extra_energy_variables(icell,j) *= internal_energy_factor;
			}
#else
			cell_gas_internal_energy(icell) = energy - kinetic_energy;
#endif /* num_extra_energy_variables > 0 */
		}
#endif /* !ADVECT_EXTRA_ENERGIES */
	}

#ifdef COOLING
#ifdef RADIATIVE_TRANSFER
	start_time( RT_COOLING_TIMER );
	rtApplyCooling(level,num_level_cells,level_cells);
	end_time( RT_COOLING_TIMER );
#else
	start_time( COOLING_TIMER );
	hydro_apply_cooling(level,num_level_cells,level_cells);
	end_time( COOLING_TIMER );
#endif /* RADIATIVE_TRANSFER */
#else
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	start_time(COOLING_TIMER);
	hydro_apply_electron_heating(level,num_level_cells,level_cells); 
	end_time(COOLING_TIMER);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* COOLING */

#ifdef EXTRA_PRESSURE_SOURCE
	hydro_zero_extra_source_vars(level,num_level_cells,level_cells);
#endif /* EXTRA_PRESSURE_SOURCE */

	cart_free( level_cells );

	end_time( WORK_TIMER );
}

void apply_hydro_fluxes( int icell, double factor, double dxi_factor, double f[num_hydro_vars-1] ) {
	int j;

	backup_hvar(icell,0) += factor*f[0];
	backup_hvar(icell,1) += factor*f[4];
	backup_hvar(icell,mj3) += factor*f[1];
	backup_hvar(icell,mj4) += factor*f[2];
	backup_hvar(icell,mj5) += factor*f[3];
	backup_hvar(icell,5) += factor*f[5];
	ref[icell] += factor*f[6]*dxi_factor;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	backup_hvar(icell,6) += factor*f[7];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_extra_energy_variables; j++ ) {
	    backup_hvar(icell,j+6+num_electronion_noneq_vars) += factor*f[7+num_electronion_noneq_vars+j];
	}
	    
	for ( j = num_hydro_vars-num_chem_species-2; j < num_hydro_vars-2; j++ ) {
		backup_hvar(icell,j) += factor*f[j+1];
	}
}


void compute_hydro_fluxes( int cell_list[4], double f[num_hydro_vars-1] ) {
    int i,j, irl;
	double v[num_state_vars][4]; /* not column-major order. */
	double c[2];

#ifdef GRAVITY_IN_RIEMANN
/* 	double g[4]; # The correct thing to do is g[4] with slope limiter in Riemann*/
	double g[2];
#endif

	int L2 = cell_list[0];
	int L1 = cell_list[1];
	int R1 = cell_list[2];
	int R2 = cell_list[3];

	cart_assert( cell_is_leaf(L1) && cell_is_leaf(R1) );

	for ( i = 0; i < 4; i++ ) {
		irl = cell_list[i];

		v[0][i] = cell_gas_density(irl);
		v[1][i] = MAX(cell_gas_pressure(irl),1e-30); /* Pressure floor is applied *after* gamma_eff calculation */
		for ( j = 0; j < num_extra_energy_variables; j++ ) {
			v[1][i] += cell_extra_energy_pressure(irl,j);
		}
		v[2][i] = cell_momentum(irl,j3)/cell_gas_density(irl);
		v[3][i] = cell_momentum(irl,j4)/cell_gas_density(irl);
		v[4][i] = cell_momentum(irl,j5)/cell_gas_density(irl);
		/* gamma_eff = (g1*P1+g2*P2+...)/(P1+P2+...) */
		v[5][i] = cell_gas_gamma(irl)*MAX(cell_gas_pressure(irl),1e-30);
		for ( j = 0; j < num_extra_energy_variables; j++ ) {
			v[5][i] += extra_energy_gamma(j)*cell_extra_energy_pressure(irl,j);
		}
		v[5][i] /= v[1][i];

#ifdef EXTRA_PRESSURE_SOURCE
		v[1][i] += cell_extra_pressure_source(irl);
#endif /* EXTRA_PRESSURE_SOURCE */
		v[1][i] = MAX( pressure_floor * v[0][i]*v[0][i], v[1][i]);
		v[6][i] = constants->gamma;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
		v[7][i] = cell_electron_internal_energy(irl);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
		for ( j = 0; j < num_extra_energy_variables; j++ ) {
			v[j+7+num_electronion_noneq_vars][i] = cell_extra_energy_variables(irl,j);
		}
		for ( j = 0; j < num_chem_species; j++ ) {
			v[num_hydro_vars-num_chem_species-1+j][i] = cell_advected_variable(irl,j)/cell_gas_density(irl);
		}
#if defined(EXTRA_PRESSURE_SOURCE) || num_extra_energy_variables > 0
		v[num_hydro_vars-1][i] = cell_gas_internal_energy(irl);
#endif /* EXTRA_PRESSURE_SOURCE || num_extra_energy_variables > 0 */

#ifdef GRAVITY_IN_RIEMANN
		/* Roughly truelove 98 (eq 34,36) */
		/* but they want s(n-1/2) for predictor then s(n+1/2) for update. */
		if(irl==1){g[0] = 0.5*cell_accel( irl, j3 ); }
		if(irl==2){g[1] = 0.5*cell_accel( irl, j3 ); }
#endif
	}

	if ( cell_level(R1) > cell_level(L1) ) {
		c[0] = 1.0/1.5;
		c[1] = 1.0/1.25;
	} else if ( cell_level(R1) < cell_level(L1) ) {
		c[0] = 1.0/1.25;
		c[1] = 1.0/1.5;
	} else {
		if ( cell_level( L2 ) == cell_level(L1) ) {
			c[0] = 1.0;
		} else {
			c[0] = 1.0/1.25;
		}

		if ( cell_level( R2 ) == cell_level(L1) ) {
			c[1] = 1.0;
		} else {
			c[1] = 1.0/1.25;
		}
	} 

#ifdef GRAVITY_IN_RIEMANN
	fluxh( dtx, dtx2, v, g, c, f );
#else	
	fluxh( dtx, dtx2, v, c, f );
#endif

#if defined(SGS_TURBULENCE) && !defined(SGST_ISOTROPIC)
	hydro_sgst_add_turbulent_fluxes( cell_list, f );
#endif /* defined(SGS_TURBULENCE) && !defined(SGST_ISOTROPIC) */

	if(apply_lapidus_viscosity) lapidus( dtx2, L1, R1, sweep_direction, j3, j4, j5, v, f );
}

void hydro_copy_vars( int level, int direction ) {
	int i, j;
	int icell;
	int num_level_cells;
	int *level_cells;

#if nDim != 3
	#error	hydro_copy_vars only works for nDim = 3
#endif

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );

	switch(direction)
	  {
	  case HYDRO_COPY_ALL:
	    {
#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,cell_child_oct,cell_vars_data,direction,backup_hvars,backup_dirty,ref,remaining_density)
	    for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			backup_hvar(icell,0) = cell_gas_density(icell);
			backup_hvar(icell,1) = cell_gas_energy(icell);
			backup_hvar(icell,2) = cell_momentum(icell,0);
			backup_hvar(icell,3) = cell_momentum(icell,1);
			backup_hvar(icell,4) = cell_momentum(icell,2);
			backup_hvar(icell,5) = cell_gas_internal_energy(icell);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			backup_hvar(icell,6) = cell_electron_internal_energy(icell);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

			for ( j = 0; j < num_extra_energy_variables; j++ ) {
			    backup_hvar(icell,j+6+num_electronion_noneq_vars) = cell_extra_energy_variables(icell,j);
			}
			for ( j = 0; j < num_chem_species; j++ ) {
			  backup_hvar(icell,num_hydro_vars-num_chem_species-2+j) = cell_advected_variable(icell,j);
			}

			ref[icell] = 0.0;
			backup_dirty[icell] = 0;
#ifdef HYDRO_TRACERS_MONTE_CARLO
			remaining_density[icell] = cell_gas_density(icell);
#endif /* HYDRO_TRACERS_MONTE_CARLO */
		}
		break;
	    }
	  case HYDRO_RESTORE_ALL:
	    {


#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,cell_child_oct,cell_vars_data,direction,backup_hvars,ref)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_gas_density(icell) = MAX( 1.0e-30, backup_hvar(icell,0) );
			cell_gas_energy(icell) = MAX( 1.0e-30, backup_hvar(icell,1) );
			cell_momentum(icell,0) = backup_hvar(icell,2);
			cell_momentum(icell,1) = backup_hvar(icell,3);
			cell_momentum(icell,2) = backup_hvar(icell,4);
			cell_gas_internal_energy(icell) = MAX( 1.0e-30, backup_hvar(icell,5) );
                                        
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			cell_electron_internal_energy(icell) = backup_hvar(icell,6);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

			for ( j = 0; j < num_extra_energy_variables; j++ ) {
		    	        cell_extra_energy_variables(icell,j) = MAX( 0.0, backup_hvar(icell,j+6+num_electronion_noneq_vars));
			}
			for ( j = 0; j < num_chem_species; j++ ) {
				cell_advected_variable(icell,j) = MAX( 1.0e-30, 
						backup_hvar(icell,num_hydro_vars-num_chem_species-2+j) );
			}

		}
		break;
	    }
	  case HYDRO_RESTORE_CLEAN:
	    {
#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,cell_child_oct,cell_vars_data,direction,backup_hvars,ref,backup_dirty,remaining_density)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( !backup_dirty[icell] ) {                                                                                                           
				cell_gas_density(icell) = MAX( 1.0e-30, backup_hvar(icell,0) );
				cell_gas_energy(icell) = MAX( 1.0e-30, backup_hvar(icell,1) );
				cell_momentum(icell,0) = backup_hvar(icell,2);
				cell_momentum(icell,1) = backup_hvar(icell,3);
				cell_momentum(icell,2) = backup_hvar(icell,4);
				cell_gas_internal_energy(icell) = MAX( 1.0e-30, backup_hvar(icell,5) );
                                        
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				cell_electron_internal_energy(icell) = backup_hvar(icell,6);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

				for ( j = 0; j < num_extra_energy_variables; j++ ) {
				        cell_extra_energy_variables(icell,j) = MAX( 0.0, backup_hvar(icell,j+6+num_electronion_noneq_vars)); 
				}
				for ( j = 0; j < num_chem_species; j++ ) {
					cell_advected_variable(icell,j) = MAX( 1.0e-30, 
							backup_hvar(icell,num_hydro_vars-num_chem_species-2+j) );
				}

#ifdef HYDRO_TRACERS_MONTE_CARLO
				remaining_density[icell] = cell_gas_density(icell);
#endif /* HYDRO_TRACERS_MONTE_CARLO */
			}
		}
		break;
	    }
	  default:
	    {
	      cart_error("Invalid hydro_copy_vars:direction parameter");
	    }
	  }		

	cart_free( level_cells );

	end_time( WORK_TIMER );
}

#endif /*HYDRO*/
