#ifndef _COSMOLOGY_INCL_
#define _COSMOLOGY_INCL_

#define Omega_M0 0.3036
#define Omega_L0 0.6964
#define Omega_b0 0.0479
#define sigma_8 0.8
#define cosmo_ns 0.968
#define cosmo_Y 0.24
#define H_z0 68.14

int init_load_cosmology_data();
double cosmo_age_z(double z);
double cosmo_z_age(double age);
double cosmo_HubbleConst_z(double z);
double cosmo_rhoc_z(double z);
double cosmo_Om_z(double z);
double cosmo_Ob_z(double z);
int free_loaded_cosmology_data();

#endif
