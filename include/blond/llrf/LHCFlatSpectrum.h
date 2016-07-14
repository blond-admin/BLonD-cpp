/*
 * LHCFlatSpectrum.h
 *
 *  Created on: June 22, 2016
 *      Author: kiliakis
 *
 * Generate LHC-type phase noise from a band-limited spectrum.
 * Input frequency band using 'fmin' and 'fmax' w.r.t. the synchrotron
 * frequency. Input double-sided spectrum amplitude [rad^2/Hz] using
 * 'initial_amplitude'. Fix seeds to obtain reproducible phase noise.
 * Select 'time_points' suitably to resolve the spectrum in frequency
 * domain. After 'corr_time' turns, the seed is changed to cut numerical
 * correlated sequences of the random number generator.
 *
 */

#ifndef LLRF_LHCFLATSPECTRUM_H_
#define LLRF_LHCFLATSPECTRUM_H_


#include <globals.h>
#include <utilities.h>
#include <configuration.h>
#include <constants.h>

class LHCFlatSpectrum {
private:

public:
   enum predistortion_t {
      exponential,
      linear,
      hyperbolic,
      weightfunction,
      None
   };

   uint fNt;
   uint fCorr;
   ftype fFMin, fFMax;
   ftype fAi;
   int fSeed1, fSeed2;
   predistortion_t fPredistortion;
   uint fNTurns;
   f_vector_t fDphi;
   f_vector_t fFs;

   LHCFlatSpectrum(uint time_points, uint corr_time = 10000,
                   ftype fmin = 0.8571, ftype fmax = 1.1,
                   ftype initial_amplitude = 1e-6,
                   int seed1 = 1234, int seed2 = 7564,
                   predistortion_t predistortion =
                      predistortion_t::None);
   ~LHCFlatSpectrum();
   void generate();
};


#endif /* LLRF_LHCFLATSPECTRUM_H_ */
