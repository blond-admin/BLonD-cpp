/*
 * math_functions.h
 *
 *  Created on: Mar 21, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_MATH_FUNCTIONS_H_
#define INCLUDES_MATH_FUNCTIONS_H_

#include <cmath>
#include "sin.h"
#include <omp.h>
#include  <cassert>
#include "utilities.h"
#include "configuration.h"
#include "fft.h"
#include <algorithm>
#include <fftw3.h>

namespace mymath {

   // Wrapper function for vdt::fast_sin
   static inline ftype fast_sin(ftype x)
   {
      return vdt::fast_sin(x);
   }


   static inline ftype fast_cos(ftype x)
   {
      return vdt::fast_sin(x + M_PI_2);
   }


// linear convolution function
   static inline void convolution(const ftype *__restrict__ signal,
                                  const uint SignalLen,
                                  const ftype *__restrict__ kernel,
                                  const uint KernelLen,
                                  ftype *__restrict__ res)
   {
      const uint size = KernelLen + SignalLen - 1;

      #pragma omp parallel for
      for (uint n = 0; n < size; ++n) {
         res[n] = 0;
         const uint kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
         const uint kmax = (n < SignalLen - 1) ? n : SignalLen - 1;
         //uint j = n - kmin;
         for (uint k = kmin; k <= kmax; k++) {
            res[n] += signal[k] * kernel[n - k];
            //--j;
         }
      }

   }


   static inline void convolution_with_ffts(f_vector_t signal,
         f_vector_t kernel,
         f_vector_t &res)
   {
      complex_vector_t v1;//(signal.size());
      complex_vector_t v2;//(kernel.size());
      const uint size = signal.size() + kernel.size() - 1;
      res.resize(size);

      fft::rfft(signal, v1, size, omp_get_max_threads());
      fft::rfft(kernel, v2, size, omp_get_max_threads());

      std::transform(v1.begin(), v1.end(),
                     v2.begin(), v1.begin(),
                     std::multiplies<complex_t>());

      fft::irfft(v1, res, size, omp_get_max_threads());

   }

   /*
   // Parameters are like python's np.interp
   // @x: x-coordinates of the interpolated values
   // @xp: The x-coords of the data points
   // @fp: the y-coords of the data points
   // @y: the interpolated values, same shape as x
   // @left: value to return for x < xp[0]
   // @right: value to return for x > xp[last]
      static inline void lin_interp(const std::vector<ftype> &x, const std::vector<ftype> &xp,
                                    const std::vector<ftype> &fp, std::vector<ftype> &y,
                                    const ftype left = 0, const ftype right = 0)
      {
         //assert(y.empty());

         gsl_interp *interp =
            gsl_interp_alloc(gsl_interp_linear, xp.size());

         gsl_interp_init(interp, &xp[0], &fp[0], xp.size());

         gsl_interp_accel *acc = gsl_interp_accel_alloc();

         y.resize(x.size());
         for (uint i = 0; i < x.size(); ++i) {
            double val;
            if (x[i] < interp->xmin) {
               //std::cout << "here\n";
               val = left;
            } else if (x[i] > interp->xmax) {
               //std::cout << "here\n";

               val = right;
            } else {
               val = gsl_interp_eval(interp, &xp[0],
                                     &fp[0], x[i],
                                     acc);
            }
            y[i] = val;
         }

         gsl_interp_free(interp);
         gsl_interp_accel_free(acc);

      }
   */

// Parameters are like python's np.interp
// @x: x-coordinates of the interpolated values
// @xp: The x-coords of the data points
// @fp: the y-coords of the data points
// @y: the interpolated values, same shape as x
// @left: value to return for x < xp[0]
// @right: value to return for x > xp[last]
   static inline void lin_interp(const std::vector<ftype> &x,
                                 const std::vector<ftype> &xp,
                                 const std::vector<ftype> &yp,
                                 std::vector<ftype> &y,
                                 const ftype left = 0.0,
                                 const ftype right = 0.0)
   {
      //assert(y.empty());
      assert(xp.size() == yp.size());

      y.resize(x.size());

      const uint N = x.size();
      //const uint M = xp.size();
      const auto max = xp.back();
      const auto min = xp.front();
      const auto end = xp.end();
      const auto begin = xp.begin();

      uint k = 0;
      while (x[k] < min and k < N) {
         y[k] = left;
         ++k;
      }

      auto j = begin + k;

      for (uint i = k; i < N; ++i) {
         if (x[i] > max) {
            y[i] = right;
            continue;
         }
         j = std::lower_bound(j, end, x[i]);
         const auto pos = j - begin;
         if (*j == x[i]) {
            y[i] = yp[pos];
         } else {
            y[i] = yp[pos - 1]
                   + (yp[pos] - yp[pos - 1])
                   * (x[i] - xp[pos - 1])
                   / (xp[pos] - xp[pos - 1]);
         }
      }

   }



// Function to implement integration of f(x) over the interval
// [a,b] using the trapezoid rule with nsub subdivisions.
   template<typename T>
   static inline ftype *cum_trapezoid(const T *f,
                                      const T deltaX,
                                      const uint nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype *psum = new ftype[nsub];
      psum[0] = 0;

      // increment the partial sum
      //#pragma omp parallel for
      for (uint i = 1; i < nsub; ++i)
         psum[i] = psum[i - 1] + (f[i] + f[i - 1]) * (deltaX / 2.0);

      return psum;

   }


   template<typename T>
   static inline ftype trapezoid(T *f,
                                 const ftype *deltaX,
                                 const uint nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions

      ftype psum = 0.0;
      // increment the partial sum
      //#pragma omp parallel for reduction(+ : psum)
      for (uint index = 1; index < nsub; ++index) {
         psum += (f[index] + f[index - 1])
                 * (deltaX[index] - deltaX[index - 1]);
      }

      // return approximation
      return psum / 2;

   }


   template<typename T>
   static inline ftype trapezoid(T *f,
                                 const ftype deltaX,
                                 const uint nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype psum = f[0] + f[nsub - 1]; //f(a)+f(b);
      //ftype deltaX = (b-a)/nsub;

      // increment the partial sum
      #pragma omp parallel for reduction(+ : psum)
      for (uint index = 1; index < nsub - 1; ++index) {
         psum += 2 * f[index];
      }

      // multiply the sum by the constant deltaX/2.0
      psum = (deltaX / 2) * psum;

      // return approximation
      return psum;

   }


   template<typename T>
   static inline uint min(T *a,
                          uint size,
                          uint step = 1)
   {
      uint p = 0;
      T min = a[0];
      //#pragma omp parallel for  shared(p) reduction(min : min)
      for (uint i = 1; i < size; i += step) {
         if (a[i] < min) {
            min = a[i];
            p = i;
         }
      }
      return p;

   }

   template<typename T>
   static inline uint max(T *a,
                          uint size,
                          uint step = 1)
   {
      uint p = 0;
      T max = a[0];
      //#pragma omp parallel for shared(p) reduction(max : max)
      for (uint i = 1; i < size; i += step) {
         if (a[i] > max) {
            max = a[i];
            p = i;
         }
      }
      return p;

   }


   static inline void linspace(ftype *a,
                               const ftype start,
                               const ftype end,
                               const uint n,
                               const uint keep_from = 0)
   {
      const ftype step = (end - start) / (n - 1);
      //ftype value = start;
      //#pragma omp parallel for
      for (uint i = 0; i < n; ++i) {
         if (i >= keep_from)
            a[i - keep_from] = start + i * step;
         //value += step;
      }
   }


   template<typename T>
   static inline std::vector<T> arange(T start,
                                       T stop,
                                       T step = 1)
   {
      std::vector<T> values;
      for (T value = start; value < stop; value += step)
         values.push_back(value);
      return values;
   }


   template<typename T>
   static inline ftype mean(const T data[],
                            const uint n)
   {
      ftype m = 0.0;
      #pragma omp parallel for reduction(+ : m)
      for (uint i = 0; i < n; ++i) {
         m += data[i];
      }
      return m / n;
   }


   template<typename T>
   static inline ftype standard_deviation(const T data[],
                                          const int n,
                                          const ftype mean)
   {
      ftype sum_deviation = 0.0;
      #pragma omp parallel for reduction(+ : sum_deviation)
      for (int i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return std::sqrt(sum_deviation / n);
   }


   template<typename T>
   static inline ftype standard_deviation(const T data[],
                                          const int n)
   {
      const ftype mean = mymath::mean(data, n);
      ftype sum_deviation = 0.0;
      #pragma omp parallel for reduction(+ : sum_deviation)
      for (int i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return std::sqrt(sum_deviation / n);
   }


   template <typename T>
   int sign(T val)
   {
      return (T(0) < val) - (val < T(0));
   }

}

#endif /* INCLUDES_MATH_FUNCTIONS_H_ */
