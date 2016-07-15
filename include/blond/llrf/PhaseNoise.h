/*
* @Author: Konstantinos Iliakis
* @Date:   2016-06-21 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-06-21 14:38:41
*/

#ifndef LLRF_PHASENOISE_H_
#define LLRF_PHASENOISE_H_

#include <blond/configuration.h>
#include <blond/utilities.h>

class API PhaseNoise {
  private:
  public:
    enum transform_t { r, c, None };

    f_vector_t fFreqArray;
    f_vector_t fRes;
    uint fResLen;
    ftype fFreqArrayMax;
    int fSeed1, fSeed2;
    uint fNt;
    ftype fDt;
    f_vector_t fT;
    f_vector_t fDphi;

    void spectrum_to_phase_noise(transform_t transform = transform_t::None);
    PhaseNoise(f_vector_t freqArray, f_vector_t realPartOfSpectrum,
               int seed1 = 0, int seed2 = 0);

    ~PhaseNoise();
};

#endif /* LLRF_PHASENOISE_H_ */
