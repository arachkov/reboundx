/**
 * General central force.
 * 
 * This example shows how to add a general central force.
 * If you have GLUT installed for the visualization, press 'w' and/or 'c' for a clearer view of
 * the whole orbit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);
double reb_E0;
double rebx_E0;
double tmax = 1.e6*2.*M_PI;

static double sqrt9(double a){
    double x = 1.;
    for (int k=0; k<20;k++){  // A smaller number should be ok too.
        double x8 = x*x*x*x*x*x*x*x;
        x += (a/x8-x)/9.;
    }
    return x;
}

int main(int argc, char* argv[]){

	int procs_N = 1;
	int procs_N_start = 0;
	if (argc>=2){
		procs_N = atoi(argv[1]);    
	}
	if (argc>=3){
		procs_N_start = atoi(argv[2]);    
	}
	int proc_i;
	for(proc_i=procs_N_start; proc_i<procs_N-1; proc_i++) {
		int pid = fork();
		if (pid == 0) {
			break;
		}
	}

	// Read initial conditions

	char filename[512];
	sprintf(filename,"notides_%04d.bin",proc_i);

    struct reb_simulation* sim = reb_create_simulation();

    struct reb_particle star = {0};
    star.m     = 1.;   
    reb_add(sim, star);
   
    double m = 1./333000;
    double a = 1.; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.0167;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    reb_add(sim, planet);

    sim->dt = 10./365.25*2.*M_PI;    
    sim->simulationarchive_interval = 2.*M_PI*5.1434563422923e4;
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->heartbeat      = heartbeat;

    reb_move_to_com(sim);

    char rebx_filename[512] = "rebx_effects_notides.bin";

    struct rebx_extras* rebx = rebx_init(sim);

    rebx_add(rebx, "moon_quadrupole_quinn");

    struct reb_particle* ps = sim->particles;
    double* f_mqq = rebx_add_param(&ps[1], "f_mqq", REBX_TYPE_DOUBLE);
    *f_mqq = 0.9473;

    double* moon_mass_mqq = rebx_add_param(&ps[1], "moon_mass_mqq", REBX_TYPE_DOUBLE);
    *moon_mass_mqq = 3.69e-8;

    double* moon_distance_mqq = rebx_add_param(&ps[1], "moon_distance_mqq", REBX_TYPE_DOUBLE);
    *moon_distance_mqq = 0.00257;
//    int* laskar_tides_mqq = rebx_add_param(&ps[1], "laskar_tides_mqq", REBX_TYPE_INT);
//    *laskar_tides_mqq = 1; // set to anything to turn tides on. comment out to set tides off.

    rebx_output_binary(rebx, rebx_filename);

    reb_E0 = reb_tools_energy(sim); // relativistic hamiltonian
    rebx_E0 = rebx_moon_quadrupole_quinn_hamiltonian(sim); // relativistic hamiltonian
	sim->simulationarchive_filename = filename;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 3.e3*2.*M_PI)){
        reb_output_timing(sim, tmax);

//        double a1 = 0.10283022841433383;
//        double a2 = -0.016869154631638014;
//        double tide_factor = sqrt9(1.0+(9.0*a1)*1.e-9*2.*M_PI*sim->t+a2*1.e-18*4.*M_PI*M_PI*sim->t*sim->t);
//        double E0 = reb_E0 + rebx_E0*tide_factor*tide_factor;
        double E0 = reb_E0 + rebx_E0;

        double Ef = rebx_moon_quadrupole_quinn_hamiltonian(sim) + reb_tools_energy(sim); // relativistic hamiltonian

        FILE* f = fopen("test_notides_relative_energy.txt","a");
        struct reb_particle com = reb_get_com(sim);
        fprintf(f,"%e %e %e %e \n",sim->t/(2.*M_PI), fabs((Ef-E0)/E0), com.x, com.y);
        fclose(f);

    }
}
