/*
 * Distributions.h
 *
 *  Created on: Mar 16, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_DISTRIBUTIONS_H_
#define BEAMS_DISTRIBUTIONS_H_

#include <cmath>
#include "globals.h"
#include "configuration.h"
#include "constants.h"
#include <stdlib.h>
#include <random>

inline void longitudinal_bigaussian(ftype sigma_dt,
                                    ftype sigma_dE = 0,
                                    int seed = 0,
                                    bool reinsertion = false)
{
   if (GP->n_sections > 1) {
      dprintf(
         "WARNING: longitudinal_bigaussian is not yet properly computed for several sections!");
   }
   if (RfP->n_rf > 1) {
      dprintf(
         "WARNING: longitudinal_bigaussian for multiple RF is not yet implemented");
   }

   uint counter = RfP->counter;
   ftype harmonic = RfP->harmonic[0][counter];
   ftype energy = GP->energy[0][counter];
   ftype beta = GP->beta[0][counter];
   ftype omega_RF = RfP->omega_RF[0][counter];
   ftype phi_s = RfP->phi_s[counter];
   ftype phi_RF = RfP->phi_RF[0][counter];

   ftype voltage, eta0, phi_b;
   if (sigma_dE == 0) {
      voltage = GP->charge * RfP->voltage[0][counter];
      eta0 = RfP->eta_0(counter);
      phi_b = omega_RF * sigma_dt + phi_s;
      sigma_dE = sqrt(
                    voltage * energy * beta * beta
                    * (cos(phi_b) - cos(phi_s)
                       + (phi_b - phi_s) * sin(phi_s))
                    / (constant::pi * harmonic * eta0));
   }

   Beam->sigma_dE = sigma_dE;
   Beam->sigma_dt = sigma_dt;

#ifdef FIXED_PARTICLES
   for (uint i = 0; i < Beam->n_macroparticles; ++i) {
      ftype r = 1.0 * (i + 1) / Beam->n_macroparticles;
      //ftype r = distribution(generator);
      Beam->dt[i] = sigma_dt * r + (phi_s - phi_RF) / omega_RF;
      //r = 1.0 * rand() / RAND_MAX;
      //r = distribution(generator);
      Beam->dE[i] = sigma_dE * r;
      //dprintf("beam_dE: %.8lf \n", Beam->dE[i]);
   }
#else
   std::default_random_engine generator(seed);
   std::normal_distribution < ftype > distribution(0.0, 1.0);
   for (uint i = 0; i < Beam->n_macroparticles; ++i) {
      //ftype r = 1.0 * rand() / RAND_MAX;
      ftype r = distribution(generator);
      Beam->dt[i] = sigma_dt * r + (phi_s - phi_RF) / omega_RF;
      //r = 1.0 * rand() / RAND_MAX;
      r = distribution(generator);
      Beam->dE[i] = sigma_dE * r;
      //dprintf("beam_dE: %.8lf \n", Beam->dE[i]);
   }
#endif

// TODO if reinsertion == true
   if (reinsertion) {
      ;
   }

}

#endif /* BEAMS_DISTRIBUTIONS_H_ */
