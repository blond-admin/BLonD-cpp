/*
* @Author: Konstantinos Iliakis
* @Date:   2016-06-22 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-06-22 14:38:41

* Methods to generate RF phase noise from noise spectrum and feedback noise
* amplitude as a function of bunch length**

* :Authors: **Helga Timko**
*/

#include "LHCFlatSpectrum.h"
#include <algorithm>
#include <constants.h>
#include <globals.h>
#include <math_functions.h>


LHCFlatSpectrum::LHCFlatSpectrum(uint time_points,
                                 uint corr_time,
                                 ftype fmin, ftype fmax,
                                 ftype initial_amplitude,
                                 int seed1, int seed2,
                                 predistortion_t predistortion)
{

   fNt = time_points;
   fCorr = corr_time;
   fFMin = fmin;
   fFMax = fmax;
   fPredistortion = predistortion;
   fSeed1 = seed1;
   fSeed2 = seed2;
   fNTurns = GP->n_turns;
   fDphi.resize(fNTurns + 1, 0);

   if (fPredistortion != predistortion_t::None) {
      // Overwrite frequencies
      fFMin = 0.8571;
      fFMax = 1.001;
   }

   if (fNt < 2 * fCorr) {
      std::cerr << "ERROR: Need more time point in LHCFlatSpectrum.\n";
      exit(-1);
   }

   // Synchrotron frequency array
   f_vector_t phis(fNTurns + 1);
   calc_phi_s(
      phis.data(),
      RfP,
      RfParameters::accelerating_systems_t::as_single);


   fFs.resize(phis.size());
   // std::cout << "c " << constant::c << "\n";
   // std::cout << "ring_circumference " << GP->ring_circumference << "\n";
   // std::cout << "harmonic[0] " << RfP->harmonic[0] << "\n";
   // std::cout << "voltage[0] " << RfP->voltage[0] << "\n";

   //uint row = RfP->section_index * (fNTurns + 1);
   for (uint i = 0; i < fFs.size(); ++i) {
      fFs[i] = constant::c
               / GP->ring_circumference
               * std::sqrt(RfP->harmonic[i] * RfP->voltage[i]
                           * std::fabs(RfP->eta_0(i) * std::cos(phis[i]))
                           / (2 * constant::pi * RfP->energy(i)));
   }

}

void LHCFlatSpectrum::generate() {}

