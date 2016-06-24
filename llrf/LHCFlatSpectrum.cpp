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
#include "PhaseNoise.h"



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
   fAi = initial_amplitude;
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
               * std::sqrt(RfP->harmonic[RfP->idx][i] * RfP->voltage[RfP->idx][i]
                           * std::fabs(RfP->eta_0(i) * std::cos(phis[i]))
                           / (2 * constant::pi * RfP->energy(i)));
   }

}



void LHCFlatSpectrum::generate()
{
   for (uint i = 0; i < fNTurns / fCorr; ++i) {

      // Scale amplitude to keep area (phase noise amplitude) constant
      auto k = i * fCorr;     // Current time step
      auto ampl = fAi * fFs[0] / fFs[k];

      // Calculate the frequency step
      uint nf = fNt / 2 + 1;      // #points in frequency domain
      auto df = GP->f_rev[k] / fNt;

      // Construct spectrum
      auto nmin = std::floor(fFMin * fFs[k] / df);
      auto nmax = std::ceil(fFMax * fFs[k] / df);

      f_vector_t freq(nf);
      mymath::linspace(freq.data(), 0, nf * df, nf);

      // std::cout << k << "\n";
      // std::cout << ampl << "\n";
      // std::cout << nf << "\n";
      // std::cout << df << "\n";
      // std::cout << nmin << "\n";
      // std::cout << nmax << "\n";

      // To compensate the noch due to PL at central frequency
      f_vector_t spectrum;

      if (fPredistortion == predistortion_t::exponential) {
         spectrum.resize(nmin, 0);

         auto v1 = mymath::arange<ftype>(0, nmax - nmin + 1);

         auto f1 = [ampl, nmax, nmin](ftype x) {
            return ampl * std::exp(std::log(100.0) * x / (nmax - nmin));
         };

         std::transform(v1.begin(),
                        v1.end(),
                        v1.begin(),
                        f1);
         // util::dump(v1, "v1\n");

         spectrum.insert(spectrum.end(), v1.begin(), v1.end());

         spectrum.resize(nf, 0);
      } else if (fPredistortion == predistortion_t::linear) {
         spectrum.resize(nmin, 0);

         f_vector_t v1(nmax - nmin + 1);
         mymath::linspace(v1.data(), 0, ampl, nmax - nmin + 1);

         spectrum.insert(spectrum.end(), v1.begin(), v1.end());

         spectrum.resize(nf, 0);
      } else if (fPredistortion == predistortion_t::hyperbolic) {
         spectrum.resize(nmin, 0);

         auto v1 = mymath::arange<ftype>(nmin, nmax + 1);

         auto f1 = [ampl, nmax, nmin](ftype x) {
            return ampl * 1.0 / (1 + 0.99 * (nmin - x) / (nmax - nmin));
         };

         std::transform(v1.begin(),
                        v1.end(),
                        v1.begin(),
                        f1);
         spectrum.insert(spectrum.end(), v1.begin(), v1.end());

         spectrum.resize(nf, 0);
      } else if (fPredistortion == predistortion_t::weightfunction) {
         spectrum.resize(nmin, 0);

         // frequency relative to fs0
         f_vector_t frel(&freq[nmin], &freq[nmax + 1]);

         std::transform(frel.begin(),
                        frel.end(),
                        frel.begin(),
                        std::bind2nd(std::divides<ftype>(), fFs[k]));

         // truncate center freqs
         auto f1 = [](ftype x) {return x > 0.999 ? 0.999 : x;};
         std::transform(frel.begin(),
                        frel.end(),
                        frel.begin(),
                        f1);
         // rms bunch length in rad corresponding to 1.2 ns
         auto sigma = 0.754;
         auto gamma = 0.577216;
         auto pi = constant::pi;
         auto f2 = [sigma, gamma, pi](ftype x) {
            auto sigma2 = sigma * sigma;
            return std::pow(4 * pi * x / sigma2, 2)
                   * std::exp(-16 * (1 - x) / sigma2)
                   + 0.25
                   * pow(1 + 8 * x / sigma2
                         * std::exp(-8 * (1 - x) / sigma2)
                         * (gamma + std::log(8 * (1 - x) / sigma2)
                            + 8 * (1 - x) / sigma2), 2);
         };

         const auto factor = ampl / f2(frel[0]);
         for (uint i = 0; i < frel.size(); ++i) {
            frel[i] = factor * f2(frel[i]);
            //frel[i] *= factor;
         }

         spectrum.insert(spectrum.end(), frel.begin(), frel.end());

         spectrum.resize(nf, 0);

      } else {
         spectrum.resize(nmin, 0);
         spectrum.resize(nmax + 1, ampl);
         spectrum.resize(nf, 0);
      }

      //util::dump(spectrum, "spectrum\n");
      auto noise = new PhaseNoise(freq, spectrum, fSeed1, fSeed2);
      noise->spectrum_to_phase_noise();
      fSeed1 += 239;
      fSeed2 += 158;

      // Fill phase noise array

      const uint kmax = i < fNTurns / fCorr - 1 ? (i + 1) * fCorr : fNTurns + 1;

      std::copy(noise->fDphi.begin(),
                noise->fDphi.begin() + kmax - k,
                fDphi.begin() + k);

      //return;

   }


}

