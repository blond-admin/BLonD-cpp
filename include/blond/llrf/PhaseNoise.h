/*
* @Author: Konstantinos Iliakis
* @Date:   2016-06-21 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-06-21 14:38:41
*/

#ifndef LLRF_PHASENOISE_H_
#define LLRF_PHASENOISE_H_

#include <blond/utilities.h>
#include <blond/configuration.h>
#include <blond/globals.h>

class API PhaseNoise {
private:

public:

   enum transform_t {
      r,
      c,
      transform_none
   };

   enum predistortion_t {
      exponential,
      linear,
      hyperbolic,
      weightfunction,
      predistortion_none
   };

   uint fCorr;
   ftype fFMin, fFMax;
   ftype fAi;
   int fSeed1, fSeed2;
   predistortion_t fPredistortion;
   uint fNTurns;
   f_vector_t fDphi;
   f_vector_t fFs;



   void spectrum_to_phase_noise(f_vector_t &t,
                                f_vector_t &dphi,
                                const f_vector_t &freq_array,
                                const f_vector_t &ReS,
                                transform_t transform =
                                   transform_t::transform_none);
   PhaseNoise() {};
   virtual ~PhaseNoise() {};
   virtual void generate() = 0;

};


class API LHCFlatSpectrum : public PhaseNoise {
private:

public:

   uint fNt;

   LHCFlatSpectrum(uint time_points, uint corr_time = 10000,
                   ftype fmin = 0.8571, ftype fmax = 1.1,
                   ftype initial_amplitude = 1e-6,
                   int seed1 = 1234, int seed2 = 7564,
                   predistortion_t predistortion =
                      predistortion_t::predistortion_none);
   ~LHCFlatSpectrum();
   void generate();

};


class API PSBPhaseNoiseInjection : public PhaseNoise {
private:

public:
   enum rescale_ampl_t {
      with_sync_freq,
      no_scaling
   };

   ftype fDeltaF;
   rescale_ampl_t fRescaleAmpl;

   PSBPhaseNoiseInjection(ftype delta_f = 1.0, uint corr_time = 10000,
                          ftype fmin = 0.8571, ftype fmax = 1.1,
                          ftype initial_amplitude = 1e-6,
                          int seed1 = 1234, int seed2 = 7564,
                          predistortion_t predistortion =
                             predistortion_t::predistortion_none,
                          rescale_ampl_t rescale_ampl =
                             rescale_ampl_t::with_sync_freq);
   ~PSBPhaseNoiseInjection();
   void generate();

};

#endif /* LLRF_PHASENOISE_H_ */
