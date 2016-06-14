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
//#include <gsl/gsl_interp.h>
//#include <gsl/gsl_errno.h>
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
   // @a: first vector
   // @b: second vector
   // @return convolution of a and b
   template<typename T>
   static inline std::vector<T> convolution(const std::vector<T> &a,
         const std::vector<T> &b)
   {
      std::vector<T> res;
      res.resize(a.size() + b.size() - 1);
      for (unsigned int i = 0; i < res.size(); i++) {
         unsigned int i1 = i;
         T temp = T();
         for (unsigned int j = 0; j < b.size(); ++j) {
            if (i1 >= 0 && i1 < a.size()) {
               temp += a[i1] * b[j];
            }
            i1--;
            res[i] = temp;
         }
      }
      return res;
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
   static inline void lin_interp(const std::vector<ftype> &x, const std::vector<ftype> &xp,
                                 const std::vector<ftype> &yp, std::vector<ftype> &y,
                                 const ftype left = 0.0, const ftype right = 0.0)
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
   static inline ftype *cum_trapezoid(const T *f, const T deltaX, const int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype *psum = new ftype[nsub];
      psum[0] = 0;

      // increment the partial sum
      //#pragma omp parallel for
      for (int i = 1; i < nsub; ++i)
         psum[i] = psum[i - 1] + (f[i] + f[i - 1]) * (deltaX / 2.0);

      return psum;

   }

   template<typename T>
   static inline ftype trapezoid(const T *f, const T *deltaX, const int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions

      ftype psum = 0.0;
      // increment the partial sum
      #pragma omp parallel for reduction(+ : psum)
      for (int index = 1; index < nsub; ++index) {
         psum += (f[index] + f[index - 1])
                 * (deltaX[index] - deltaX[index - 1]);
      }

      // return approximation
      return psum / 2;

   }

   template<typename T>
   static inline ftype trapezoid(const T *f, const T deltaX, const int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype psum = f[0] + f[nsub - 1]; //f(a)+f(b);
      //ftype deltaX = (b-a)/nsub;

      // increment the partial sum
      #pragma omp parallel for reduction(+ : psum)
      for (int index = 1; index < nsub - 1; ++index) {
         psum += 2 * f[index];
      }

      // multiply the sum by the constant deltaX/2.0
      psum = (deltaX / 2) * psum;

      // return approximation
      return psum;

   }
   template<typename T>
   static inline int min(T *a, int size, int step = 1)
   {
      int p = 0;
      T min = a[0];
      //#pragma omp parallel for  shared(p) reduction(min : min)
      for (int i = 1; i < size; i += step) {
         if (a[i] < min) {
            min = a[i];
            p = i;
         }
      }
      return p;

   }

   template<typename T>
   static inline int max(T *a, int size, int step = 1)
   {
      int p = 0;
      T max = a[0];
      //#pragma omp parallel for shared(p) reduction(max : max)
      for (int i = 1; i < size; i += step) {
         if (a[i] > max) {
            max = a[i];
            p = i;
         }
      }
      return p;

   }

   static inline void linspace(ftype *a, const ftype start, const ftype end,
                               const int n, const int keep_from = 0)
   {
      const ftype step = (end - start) / (n - 1);
      //ftype value = start;
      #pragma omp parallel for
      for (int i = 0; i < n; ++i) {
         if (i >= keep_from)
            a[i - keep_from] = start + i * step;
         //value += step;
      }
   }

   template<typename T>
   static inline ftype mean(const T data[], const int n)
   {
      ftype m = 0.0;
      #pragma omp parallel for reduction(+ : m)
      for (int i = 0; i < n; ++i) {
         m += data[i];
      }
      return m / n;
   }

   template<typename T>
   static inline ftype standard_deviation(const T data[], const int n,
                                          const ftype mean)
   {
      ftype sum_deviation = 0.0;
      #pragma omp parallel for reduction(+ : sum_deviation)
      for (int i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return std::sqrt(sum_deviation / n);
   }

   template<typename T>
   static inline ftype standard_deviation(const T data[], const int n)
   {
      const ftype mean = mymath::mean(data, n);
      ftype sum_deviation = 0.0;
      #pragma omp parallel for reduction(+ : sum_deviation)
      for (int i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return std::sqrt(sum_deviation / n);
   }

}

#endif /* INCLUDES_MATH_FUNCTIONS_H_ */
