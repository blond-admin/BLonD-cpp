/*
 * math_functions.h
 *
 *  Created on: Mar 21, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_MATH_FUNCTIONS_H_
#define INCLUDES_MATH_FUNCTIONS_H_

#include <cmath>
#include "configuration.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_errno.h>

namespace mymath {
   // Parameters are like python's numpy.fft.rfft
   // @in:  input data
   // @n:   number of points to use. If n < in.size() then the input is cropped
   //       if n > in.size() then input is padded with zeros
   // @out: the transformed array

   static inline void rfft(const std::vector<ftype> &in, const unsigned int n,
                           std::vector<complex_t> &out)
   {
      std::vector<ftype> v(in);
      v.resize(n, 0);

      gsl_fft_real_wavetable *real;
      gsl_fft_real_workspace *work;

      work = gsl_fft_real_workspace_alloc(n);
      real = gsl_fft_real_wavetable_alloc(n);

      gsl_fft_real_transform(v.data(), 1, n, real, work);
      //printf("result data %lu \n",v.size() );
      // unpack result into complex format
      out.clear();
      out.push_back(complex_t(v[0], 0));
      for (unsigned int i = 1; i < v.size(); i += 2) {
         out.push_back(complex_t(v[i], v[i + 1]));
      }
      //out.push_back(complex_t(v.back(), 0));

      gsl_fft_real_wavetable_free(real);
      gsl_fft_real_workspace_free(work);
   }

   // Parameters are like python's numpy.fft.irfft
   // @in:  input data
   // @n:   number of points to use. If n < in.size() then the input is cropped
   //       if n > in.size() then input is padded with zeros
   // @out: the inverse Fourier transform of input data
   // out size is 2*(m-1), where m is input size

   static inline void irfft(const std::vector<ftype> &in, const unsigned int n,
                            std::vector<complex_t> &out)
   {
      std::vector<ftype> v(in);
      v.resize(n, 0);

      gsl_fft_halfcomplex_wavetable *hc;
      gsl_fft_real_workspace *work;

      work = gsl_fft_real_workspace_alloc(n);
      hc = gsl_fft_halfcomplex_wavetable_alloc(n);

      gsl_fft_halfcomplex_inverse(v.data(), 1, n, hc, work);
      //printf("result data %lu \n",v.size() );
      // unpack result into complex format
      out.clear();
      out.push_back(complex_t(v[0], 0));
      for (unsigned int i = 1; i < v.size(); i += 2) {
         out.push_back(complex_t(v[i], v[i + 1]));
      }
      //out.push_back(complex_t(v.back(), 0));

      gsl_fft_halfcomplex_wavetable_free(hc);
      gsl_fft_real_workspace_free(work);
   }



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
      assert(y.empty());

      gsl_interp *interp =
         gsl_interp_alloc(gsl_interp_linear, xp.size());

      gsl_interp_init(interp, &xp[0], &fp[0], xp.size());

      gsl_interp_accel *acc = gsl_interp_accel_alloc();

      for (uint i = 0; i < x.size(); ++i) {
         double val;
         if (x[i] < interp->xmin) {
            val = left;
         } else if (x[i] > interp->xmax) {
            val = right;
         } else {
            val = gsl_interp_eval(interp, &xp[0],
                                  &fp[0], x[i],
                                  acc);
         }
         y.push_back(val);
      }

      gsl_interp_free(interp);
      gsl_interp_accel_free(acc);



   }


// Function to implement integration of f(x) over the interval
// [a,b] using the trapezoid rule with nsub subdivisions.
   static inline ftype *cum_trapezoid(ftype *f, ftype deltaX, int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype *psum = new ftype[nsub];
      psum[0] = 0;

      // increment the partial sum
      for (int i = 1; i < nsub; i++)
         psum[i] = psum[i - 1] + (f[i] + f[i - 1]) * (deltaX / 2.0);

      return psum;

   }

   template<typename T>
   static inline ftype trapezoid(T *f, ftype *deltaX, int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions

      ftype psum = 0;
      // increment the partial sum
      for (int index = 1; index < nsub; index++) {
         psum = psum
                + (f[index] + f[index - 1])
                * (deltaX[index] - deltaX[index - 1]);
      }

      // return approximation
      return psum / 2;

   }

   template<typename T>
   static inline ftype trapezoid(T *f, ftype deltaX, int nsub)
   {
      // initialize the partial sum to be f(a)+f(b) and
      // deltaX to be the step size using nsub subdivisions
      ftype psum = f[0] + f[nsub - 1]; //f(a)+f(b);
      //ftype deltaX = (b-a)/nsub;

      // increment the partial sum
      for (int index = 1; index < nsub - 1; index++) {
         psum = psum + 2 * f[index];
      }

      // multiply the sum by the constant deltaX/2.0
      psum = (deltaX / 2) * psum;

      // return approximation
      return psum;

   }
   template<typename T>
   static inline int min(T *a, int size, int step)
   {
      int p = 0;
      T min = a[0];
      for (int i = 0; i < size; i += step) {
         if (a[i] < min) {
            min = a[i];
            p = i;
         }
      }
      return p;

   }

   template<typename T>
   static inline int max(T *a, int size, int step)
   {
      int p = 0;
      ftype max = a[0];
      for (int i = 0; i < size; i += step) {
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
      ftype step = (end - start) / (n - 1);
      ftype value = start;
      for (int i = 0; i < n; ++i) {
         if (i >= keep_from)
            a[i - keep_from] = value;
         value += step;
      }
   }

   static inline ftype mean(const ftype data[], const int n)
   {
      ftype m = 0;
      for (int i = 0; i < n; ++i) {
         m += data[i];
      }
      return m / n;
   }

   static inline ftype standard_deviation(const ftype data[], const int n,
                                          const ftype mean)
   {
      ftype sum_deviation = 0.0;
      int i;
      for (i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return sqrt(sum_deviation / n);
   }

   static inline ftype standard_deviation(const ftype data[], const int n)
   {
      ftype mean = mymath::mean(data, n);
      ftype sum_deviation = 0.0;
      int i;
      for (i = 0; i < n; ++i)
         sum_deviation += (data[i] - mean) * (data[i] - mean);
      return sqrt(sum_deviation / n);
   }

}

#endif /* INCLUDES_MATH_FUNCTIONS_H_ */
