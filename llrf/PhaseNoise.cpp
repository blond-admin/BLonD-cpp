/*
* @Author: Konstantinos Iliakis
* @Date:   2016-06-21 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-06-21 14:38:41

* Methods to generate RF phase noise from noise spectrum and feedback noise
* amplitude as a function of bunch length**

* :Authors: **Helga Timko**
*/

#include "PhaseNoise.h"
#include <random>
#include <algorithm>
#include "constants.h"
#include "fft.h"
#include "math_functions.h"


PhaseNoise::PhaseNoise(f_vector_t freqArray,
                       f_vector_t realPartOfSpectrum,
                       int seed1, int seed2)
{

   fFreqArray = freqArray;
   fRes = realPartOfSpectrum;
   fSeed1 = seed1;
   fSeed2 = seed2;
   fResLen = fRes.size();
   fFreqArrayMax = fFreqArray.back();
   fNt = 0;
   fDt = 0;
}

// TODO test this function
void PhaseNoise::spectrum_to_phase_noise(PhaseNoise::transform_t transform)
{

   // Resolution in time domain

   if (transform == transform_t::None or transform == transform_t::r) {
      fNt = 2 * (fResLen - 1);
      fDt = 1 / (2 * fFreqArrayMax);
   } else if (transform == transform_t::c) {
      fNt = fResLen;
      fDt = 1 / fFreqArrayMax;
   } else {
      std::cerr << "ERROR: The choise of Fourier transform for the\n"
                << "RF noise generation could not be recognized.\n"
                << "Use 'r' or 'c'.\n";
   }


   // Generate white noise in time domain
   std::default_random_engine gen(fSeed1);
   std::uniform_real_distribution<> dist(0, 1);

   f_vector_t r1(fNt);
   for (auto &v : r1)
      v = dist(gen);

   gen.seed(fSeed2);

   f_vector_t r2(fNt);
   for (auto &v : r2)
      v = dist(gen);

   // for (uint i = 1; i < fNt + 1; ++i) {
   //    r1[i - 1] = r2[i - 1] = (1.0 * i) / (fNt + 1);
   // }

   //std::cout << "mean r1 : " << mymath::mean(r1.data(), r1.size()) << "\n";
   //std::cout << "mean r2 : " << mymath::mean(r2.data(), r2.size()) << "\n";

   if (transform == transform_t::None or transform == transform_t::r) {
      f_vector_t Gt(fNt);

      auto factor = 2 * constant::pi;
      auto f1 = [factor](ftype x) {return std::cos(factor * x);};
      std::transform(r1.begin(),
                     r1.end(),
                     r1.begin(),
                     f1);

      auto f2 = [](ftype x) {return std::sqrt(-2 * std::log(x));};
      std::transform(r2.begin(),
                     r2.end(),
                     r2.begin(),
                     f2);


      std::transform(r1.begin(),
                     r1.end(),
                     r2.begin(),
                     Gt.begin(),
                     std::multiplies<ftype>());

      //std::cout << "mean Gt : " << mymath::mean(Gt.data(), Gt.size()) << "\n";

      // FFT to frequency domain
      complex_vector_t Gf;
      fft::rfft(Gt, Gf);

      // auto sum = 0.0;
      // for (const auto &v : Gf)
      //    sum += std::abs(v);

      //std::cout << "mean abs(Gf) : " << sum / Gf.size() << "\n";

      // Multiply by desired noise probability density
      factor = 2 * fFreqArrayMax;
      auto f3 = [factor](ftype x) {return std::sqrt(factor * x);};

      r1.resize(fRes.size());
      std::transform(fRes.begin(),
                     fRes.end(),
                     r1.begin(),
                     f3);

      //std::cout << "mean s : " << mymath::mean(r1.data(), r1.size()) << "\n";


      auto f4 = [](ftype r, complex_t z) {return r * z;};
      std::transform(r1.begin(),
                     r1.end(),
                     Gf.begin(),
                     Gf.begin(),
                     f4);
      // sum = 0.0;
      // for (const auto &v : Gf)
      //    sum += std::abs(v);

      // std::cout << "mean abs(dPf) : " << sum / Gf.size() << "\n";

      // fft back to time domain to get final phase shift
      Gt.clear();
      fft::irfft(Gf, Gt);

      //std::cout << "mean dpt : " << mymath::mean(Gt.data(), Gt.size()) << "\n";

      // Use only real part of the phase shift and normalize
      fT.resize(fNt);
      mymath::linspace(fT.data(), 0, fNt * fDt, fNt);

      auto f6 = [](complex_t z) {return z.real();};

      fDphi.resize(Gt.size());

      std::transform(Gt.begin(),
                     Gt.end(),
                     fDphi.begin(),
                     f6);

   } else if (transform == transform_t::c) {

      complex_vector_t Gt(fNt);

      const complex_t factor = 2 * constant::pi * complex_t(0, 1);
      auto f1 = [factor](ftype x) {return std::exp(factor * x);};
      std::transform(r1.begin(),
                     r1.end(),
                     Gt.begin(),
                     f1);

      auto f2 = [](ftype x) {return std::sqrt(-2 * std::log(x));};
      std::transform(r2.begin(),
                     r2.end(),
                     r2.begin(),
                     f2);

      auto f3 = [](complex_t a, ftype b) {return a * b;};
      std::transform(Gt.begin(),
                     Gt.end(),
                     r2.begin(),
                     Gt.begin(),
                     f3);

      // FFT to frequency domain
      complex_vector_t Gf;
      fft::fft(Gt, Gf);

      // Multiply by desired noise probability density

      auto factor2 = fFreqArrayMax;
      auto f4 = [factor2](ftype x) {return std::sqrt(factor2 * x);};

      r1.resize(fRes.size());
      std::transform(fRes.begin(),
                     fRes.end(),
                     r1.begin(),
                     f4);

      auto f5 = [](ftype r, complex_t z) {return r * z;};
      std::transform(r1.begin(),
                     r1.end(),
                     Gf.begin(),
                     Gf.begin(),
                     f5);

      // fft back to time domain to get final phase shift
      Gt.clear();
      fft::ifft(Gf, Gt);

      // Use only real part of the phase shift and normalize
      fT.resize(fNt);
      mymath::linspace(fT.data(), 0, fNt * fDt, fNt);

      auto f6 = [](complex_t z) {return z.real();};

      fDphi.resize(Gt.size());

      std::transform(Gt.begin(),
                     Gt.end(),
                     fDphi.begin(),
                     f6);
   }



}