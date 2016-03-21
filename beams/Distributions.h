/*
 * Distributions.h
 *
 *  Created on: Mar 16, 2016
 *      Author: kiliakis
 */



#ifndef BEAMS_DISTRIBUTIONS_H_
#define BEAMS_DISTRIBUTIONS_H_


#include <cmath>
#include "Beams.h"
#include "../input_parameters/GeneralParameters.h"
#include "../input_parameters/RfParameters.h"
#include "../includes/configuration.h"
#include "../includes/constants.h"
#include <stdlib.h>
#include <random>

inline void longitudinal_bigaussian(GeneralParameters *gp,
                                    RfParameters *rf_params, Beams *beam, ftype sigma_dt,
                                    ftype sigma_dE = 0, int seed = 0, bool reinsertion = false) {
	if (gp->n_sections > 1) {
		dprintf(
		    "WARNING: longitudinal_bigaussian is not yet properly computed for several sections!");
	}
	if (rf_params->n_rf > 1) {
		dprintf(
		    "WARNING: longitudinal_bigaussian for multiple RF is not yet implemented");
	}

	int counter = rf_params->counter;
	ftype harmonic = rf_params->harmonic[counter];
	ftype energy = gp->energy[counter];
	ftype beta = gp->beta[counter];
	ftype omega_RF = rf_params->omega_RF[counter];
	ftype phi_s = rf_params->phi_s[counter];
	ftype phi_RF = rf_params->phi_RF[counter];

	ftype voltage, eta0, phi_b;
	if (sigma_dE == 0) {
		voltage = gp->charge * rf_params->voltage[0 + counter];
		eta0 = rf_params->eta_0(counter);
		phi_b = omega_RF * sigma_dt + phi_s;
		sigma_dE = sqrt(
		               voltage * energy * beta * beta
		               * (cos(phi_b) - cos(phi_s)
		                  + (phi_b - phi_s) * sin(phi_s))
		               / (pi * harmonic * eta0));
		//dprintf("cos(phi_s): %.12lf \n", cos(phi_s));
		//dprintf("cos(phi_b): %.12lf \n", cos(phi_b));
		//dprintf("sin(phi_s): %.12lf \n", sin(phi_s));


	}
	//dprintf("omega_RF: %.12lf \n", omega_RF);
	//dprintf("sigma_dt: %.12lf \n", sigma_dt);
	//dprintf("phi_s: %.12lf \n", phi_s);
	
	//dprintf("phi_b: %.12lf \n", phi_b);
	//dprintf("sigma_dE: %.12lf \n", sigma_dE);

	beam->sigma_dE = sigma_dE;
	beam->sigma_dt = sigma_dt;
	//srand(seed);
	#ifdef FIXED_PARTICLES
	for (int i = 0; i < beam->n_macroparticles; ++i) {
		ftype r = 1.0 * (i+1)/beam->n_macroparticles;
		//ftype r = distribution(generator);
		beam->dt[i] = sigma_dt * r + (phi_s - phi_RF) / omega_RF;
		//r = 1.0 * rand() / RAND_MAX;
		//r = distribution(generator);
		beam->dE[i] = sigma_dE * r;
		//dprintf("beam_dE: %.8lf \n", beam->dE[i]);

	}


	#else

	std::default_random_engine generator(seed);
  	std::normal_distribution<ftype> distribution(0.0,1.0);
	for (int i = 0; i < beam->n_macroparticles; ++i) {
		//ftype r = 1.0 * rand() / RAND_MAX;
		ftype r = distribution(generator);
		beam->dt[i] = sigma_dt * r + (phi_s - phi_RF) / omega_RF;
		//r = 1.0 * rand() / RAND_MAX;
		r = distribution(generator);
		beam->dE[i] = sigma_dE * r;
		//dprintf("beam_dE: %.8lf \n", beam->dE[i]);

	}

	#endif

	// TODO if reinsertion == true
	if (reinsertion)
	{
		;
	}

}

#endif /* BEAMS_DISTRIBUTIONS_H_ */
