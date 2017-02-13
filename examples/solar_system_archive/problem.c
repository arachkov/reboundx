/**
 * Simulation Archive
 *
 * This example shows how to use the Simulation Archive.
 * We integrate a two planet system forward in time using
 * the WHFast integrator. The simulation can be interrupted
 * at any time. On the next run, the program will try to reload
 * the latest data from the Simulation Archive. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "rebound.h"
#include "reboundx.h"
#include <time.h>


void heartbeat(struct reb_simulation* r);
int main(int argc, char* argv[]) {

	double cm_to_au = 6.68459e-14;

	srand ( time(NULL) );
	double random_number=((double)rand()/(double)RAND_MAX);

	double kappa = 0.75 + random_number*0.45;

	double k = 1.0;

	char filename[512];
	char rebx_filename[512];

	time_t t = time(NULL);
	struct tm tm = *localtime(&t);

	sprintf(filename,"/data_local/arachkov/solar_system/kappas/kappa_%e_date_%d_%d_%d_time_%dh_%dmin_%dsec.bin",kappa,tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
	sprintf(rebx_filename,"/data_local/arachkov/solar_system/kappas/rebx_date_%d_%d_%d_time_%dh_%dmin_%dsec.bin",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);

	// Trying to restart from the Simulation Archive.
	struct reb_simulation* r = reb_create_simulation_from_simulationarchive(filename);

	// Check if that was successful
	if (r==NULL){
		printf("No simulation archive found. Creating new simulation.\n");
		r= reb_create_simulation();
		struct reb_particle sun = {.m=1.};
		reb_add(r, sun);
		struct reb_particle mercury = reb_tools_orbit2d_to_particle(r->G, sun, 1./6023600, 3.870982252717257E-01+3.8*k*cm_to_au, kappa*2.056302512089075E-01, 2.912428058698772E+01*M_PI/180., 1.751155303115542E+02*M_PI/180.);
		reb_add(r, mercury);
		struct reb_particle venus = reb_tools_orbit2d_to_particle(r->G, sun, 1./408523.72, 7.233268496749391E-01, kappa*6.755697267164094E-03, 5.518541455452200E+01*M_PI/180., 4.990452231866427E+01*M_PI/180.);
		reb_add(r, venus);
		struct reb_particle earth = reb_tools_orbit2d_to_particle(r->G, sun, 3.003297890315729e-06, 1.000371833989169E+00, kappa*1.704239716781501E-02, 2.977668064579176E+02*M_PI/180., 3.581260865454548E+02*M_PI/180.);
		reb_add(r, earth);
		struct reb_particle mars = reb_tools_orbit2d_to_particle(r->G, sun, 1./3098708, 1.523678184302188E+00, kappa*9.331460653723893E-02, 2.865373577554387E+02*M_PI/180., 2.302024685336155E+01*M_PI/180.);
		reb_add(r, mars);
		struct reb_particle jupiter = reb_tools_orbit2d_to_particle(r->G, sun, 0.0009545325625181037, 5.205108604506466E+00, kappa*4.892306471604416E-02, 2.751196839758603E+02*M_PI/180., 2.063463646284857E+01*M_PI/180.);
		reb_add(r, jupiter);
		struct reb_particle saturn = reb_tools_orbit2d_to_particle(r->G, sun, 0.00028579654259598984, 9.581451990386764E+00, kappa*5.559928883801366E-02, 3.359006493683225E+02*M_PI/180., 3.160917714975463E+02*M_PI/180.);
		reb_add(r, saturn);
		struct reb_particle uranus = reb_tools_orbit2d_to_particle(r->G, sun, 4.365520702584404e-05, 1.922994520785785E+01, kappa*4.439340361752947E-02, 9.661124460893427E+01*M_PI/180., 1.458440916327605E+02*M_PI/180.);
		reb_add(r, uranus);
		struct reb_particle neptune = reb_tools_orbit2d_to_particle(r->G, sun, 5.149999195391201e-05, 3.009697072395906E+01, kappa*1.114818186443456E-02, 2.668275286227091E+02*M_PI/180., 2.653252378278363E+02*M_PI/180.);
		reb_add(r, neptune);
		struct reb_particle pluto = reb_tools_orbit2d_to_particle(r->G, sun, 6.572648128479933e-09, 3.950092123894740E+01, kappa*2.478618527514649E-01, 1.151532923291780E+02*M_PI/180., 2.385746094416748E+01*M_PI/180.);
		reb_add(r, pluto);
		reb_move_to_com(r);
		r->dt = 8./365.25*2.*M_PI;                      // 8 days, where G=1
		r->simulationarchive_interval = 2.*M_PI*1.e3; // output data every 100000 years
		r->ri_whfast.safe_mode = 0;                      
		r->ri_whfast.corrector = 5;    
		r->integrator = REB_INTEGRATOR_WHFAST;
		r->heartbeat  = heartbeat;
	}else{
		printf("Found simulation archive. Loaded snapshot at t=%.16f.\n",r->t);
	}

	struct rebx_extras* rebx = rebx_create_extras_from_binary(r, rebx_filename);
	if (rebx==NULL){
		printf("No reboundx simulation archive found. Creating new rebx effects.\n");
		rebx = rebx_init(r);
		struct rebx_effect* gr_params = rebx_add(rebx, "gr_potential");
		double* c = rebx_add_param(gr_params, "c", REBX_TYPE_DOUBLE);
		*c = REBX_C;
		rebx_output_binary(rebx, rebx_filename);
	}else{
		printf("Found rebx effects archive.");
	}

	r->simulationarchive_filename = filename;
	reb_integrate(r, 2.*M_PI*5e9); //time
	rebx_free(rebx);
	reb_free_simulation(r);

}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 1.e3*2.*M_PI)){
        reb_output_timing(sim, 2.*M_PI*5e9);
    }
}
