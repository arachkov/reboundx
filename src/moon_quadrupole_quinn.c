/** * @file moon_quadrupole_quinn.c
 * @brief   Models Earth's moon as a quadrupole around the Earth interacting with the Sun
 * @author  Aleksandar Rachkov, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Aleksandar Rachkov, Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Multipoles$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 A. Rachkov, D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                Quinn et al., 1991, Laskar & Gastineau, 2009 
 * C Example               :ref:`c_example_moon_quadrupole_quinn`
 * Python Example          `MoonQuadrupoleQuinn.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/MoonQuadrupoleQuinn.ipynb>`_.
 * ======================= ===============================================
 * 
 * This effect adds a correction term for
 * Add description here of effect 
 *
 * **Effect Parameters**
 * 
 * None  
 *
 * **Particle Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * m_ratio_planetmoon_mql (double)             Yes         Earth mass over Moon mass ratio.
 * a0_mql (double)         Yes         a0 parameter from the tidal dissipation model, equation 1 in Supplementary Material from Laskar & Gastineau, 2009.
 * a1_mql (double)         Yes         a1 parameter from the tidal dissipation model, equation 1 in Supplementary Material from Laskar & Gastineau, 2009.
 * a2_mql (double)         Yes         a2 parameter from the tidal dissipation model, equation 1 in Supplementary Material from Laskar & Gastineau, 2009.
 * alpha_mql (double)         Yes         alpha parameter from the tidal dissipation model, equation 1 in Supplementary Material Laskar & Gastineau, 2009.
 * f_mql (double)         Yes         correction factor f added due to representing the average lunar orbit as a circular ring, coming from equation 2 in Quinn et al., 1991.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

// Machine independent implementation of pow(*,9). 
static double sqrt9(double a){
    double x = 1.;
    for (int k=0; k<20;k++){  // A smaller number should be ok too.
        double x8 = x*x*x*x*x*x*x*x;
        x += (a/x8-x)/9.;
    }
    return x;
}

static void rebx_calculate_mqq_force(struct reb_simulation* const sim, const double f_mqq, const double moon_mass_mqq, double moon_distance_mqq, const int i){

    struct reb_particle* const particles = sim->particles;
    const struct reb_particle source = particles[0];
    const struct reb_particle p = particles[i];

    double a1 = 0.10283022841433383;
    double a2 = -0.016869154631638014;

    const double dx = p.x - source.x;
    const double dy = p.y - source.y;
    const double dz = p.z - source.z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    const int* const laskar_tides_mqq = rebx_get_param_check(&particles[i], "laskar_tides_mqq", REBX_TYPE_INT);

    if (laskar_tides_mqq){
        double a = 1.0+(9.0*a1)*1.e-9*2.*M_PI*sim->t+a2*1.e-18*4.*M_PI*M_PI*sim->t*sim->t;
        moon_distance_mqq *= sqrt9(a);
    }

    double massratio = p.m*moon_mass_mqq/((p.m + moon_mass_mqq)*(p.m + moon_mass_mqq));

    const double A = 1.0;
//    const double A = (3.0/4.0)*sim->G*source.m*(f_mqq)*moon_distance_mqq*moon_distance_mqq*massratio;
    const double gamma = -4.0;
    const double prefac = A*pow(r2, (gamma-1.)/2.);

    particles[i].ax += prefac*dx;
    particles[i].ay += prefac*dy;
    particles[i].az += prefac*dz;
    particles[0].ax -= p.m/source.m*prefac*dx;
    particles[0].ay -= p.m/source.m*prefac*dy;
    particles[0].az -= p.m/source.m*prefac*dz;
}

void rebx_moon_quadrupole_quinn(struct reb_simulation* const sim, struct rebx_effect* const effect){ 
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    for (int i=1; i<N_real; i++){
        const double* const f_mqq = rebx_get_param_check(&particles[i], "f_mqq", REBX_TYPE_DOUBLE);
        if (f_mqq != NULL){
            const double* const moon_mass_mqq = rebx_get_param_check(&particles[i], "moon_mass_mqq", REBX_TYPE_DOUBLE);
            if (moon_mass_mqq != NULL){
		        const double* const moon_distance_mqq = rebx_get_param_check(&particles[i], "moon_distance_mqq", REBX_TYPE_DOUBLE);
		        if (moon_distance_mqq != NULL){
                    rebx_calculate_mqq_force(sim, *f_mqq, *moon_mass_mqq, *moon_distance_mqq, i); // only calculates force if all parameters set
                }
            }
        }
    }
}

static double rebx_calculate_hamiltonian(struct reb_simulation* const sim, const double f_mqq, const double moon_mass_mqq, double moon_distance_mqq, const int i){

    double a1 = 0.10283022841433383;
    double a2 = -0.016869154631638014;
    const struct reb_particle* const particles = sim->particles;

    const struct reb_particle source = particles[0];
    const struct reb_particle p = particles[i];
    const double dx = p.x - source.x;
    const double dy = p.y - source.y;
    const double dz = p.z - source.z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    const int* const laskar_tides_mqq = rebx_get_param_check(&particles[i], "laskar_tides_mqq", REBX_TYPE_INT);

    if (laskar_tides_mqq){
        double a = 1.0+(9.0*a1)*1.e-9*2.*M_PI*sim->t+a2*1.e-18*4.*M_PI*M_PI*sim->t*sim->t;
        moon_distance_mqq *= sqrt9(a);
    }

    double massratio = p.m*moon_mass_mqq/((p.m + moon_mass_mqq)*(p.m + moon_mass_mqq));

//    const double A = (-1.0/4.0)*sim->G*source.m*(f_mqq)*moon_distance_mqq*moon_distance_mqq*massratio;
    const double A = 1.0;
    const double gamma = -4.0;

//    double H = p.m*A/(r2*sqrt(r2));
    double H = -p.m*A*pow(r2, (gamma+1.)/2.)/(gamma+1.);
    return H;
}

double rebx_moon_quadrupole_quinn_hamiltonian(struct reb_simulation* const sim){ 
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=1; i<N_real; i++){
        const double* const f_mqq = rebx_get_param_check(&particles[i], "f_mqq", REBX_TYPE_DOUBLE);
        if (f_mqq != NULL){
            const double* const moon_mass_mqq = rebx_get_param_check(&particles[i], "moon_mass_mqq", REBX_TYPE_DOUBLE);
            if (moon_mass_mqq != NULL){
		        const double* const moon_distance_mqq = rebx_get_param_check(&particles[i], "moon_distance_mqq", REBX_TYPE_DOUBLE);
		        if (moon_distance_mqq != NULL){
                    Htot += rebx_calculate_hamiltonian(sim, *f_mqq, *moon_mass_mqq, *moon_distance_mqq, i);
                }
            }
        }
    }
    return Htot;
}
