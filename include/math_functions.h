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
#include  <cassert>
#include "configuration.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_fft_real.h>
//#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_errno.h>

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

   static inline void real_to_complex(const std::vector<ftype> &in,
                                      std::vector<complex_t> &out)
   {
      assert(out.empty());
      for (unsigned int i = 0; i < in.size(); i++) {
         //complex_t z(in[i], 0);
         out.push_back(complex_t(in[i], 0));
      }

   }

   static inline void pack_to_complex(const std::vector<ftype> &in,
                                      std::vector<complex_t> &out)
   {
      assert(out.empty());
      for (unsigned int i = 0; i < in.size(); i += 2) {
         complex_t z(in[i], in[i + 1]);
         out.push_back(complex_t(in[i], in[i + 1]));
      }

   }


   static inline void complex_to_real(const std::vector<complex_t> &in,
                                      std::vector<ftype> &out)
   {
      assert(out.empty());

      for (complex_t z : in) {
         out.push_back(z.real());
         //out.push_back(z.imag());
      }

   }

   static inline void unpack_complex(const std::vector<complex_t> &in,
                                     std::vector<ftype> &out)
   {
      //assert(out.empty());
      for (unsigned int i = 0; i < in.size(); ++i) {
         out[2 * i] = in[i].real();
         out[2 * i + 1] = in[i].imag();
      }

   }


   // Parameters are like python's numpy.fft.rfft
   // @in:  input data
   // @n:   number of points to use. If n < in.size() then the input is cropped
   //       if n > in.size() then input is padded with zeros
   // @out: the transformed array

   static inline void rfft(std::vector<ftype> &in, const unsigned int n,
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
      //v.back() = 0;
      for (unsigned int i = 1; i < v.size(); i += 2) {
         out.push_back(complex_t(v[i], v[i + 1]));
      }
      //out.push_back(complex_t(v.back(), 0));
      // TODO if n is even then last element is zero ??
      if (n % 2 == 0)
         out.back() = complex_t(out.back().real(), 0);
      //printf("out size %lu \n", out.size());
      //out.push_back(complex_t(v.back(), 0));

      gsl_fft_real_wavetable_free(real);
      gsl_fft_real_workspace_free(work);
   }



   // Parameters are like python's numpy.fft.fft
   // @in:  input data
   // @n:   number of points to use. If n < in.size() then the input is cropped
   //       if n > in.size() then input is padded with zeros
   // @out: the transformed array

   static inline void fft(std::vector<complex_t> &in, const unsigned int n,
                          std::vector<complex_t> &out)
   {
      std::vector<ftype> v;
      v.resize(2 * n, 0);
      unpack_complex(in, v);

      gsl_fft_complex_wavetable *wave;
      gsl_fft_complex_workspace *work;

      wave = gsl_fft_complex_wavetable_alloc(n);
      work = gsl_fft_complex_workspace_alloc(n);

      gsl_fft_complex_forward(v.data(), 1, n, wave, work);

      //printf("ok inside\n");

      out.clear();

      pack_to_complex(v, out);

      out.resize(n, 0);

      gsl_fft_complex_wavetable_free(wave);
      //printf("ok here7\n");

      gsl_fft_complex_workspace_free(work);
      //printf("ok here8\n");

   }



   // Parameters are like python's numpy.fft.ifft
   // @in:  input data
   // @n:   number of points to use. If n < in.size() then the input is cropped
   //       if n > in.size() then input is padded with zeros
   // @out: the inverse Fourier transform of input data

   static inline void ifft(std::vector<complex_t> &in, const unsigned int n,
                           std::vector<complex_t> &out)
   {
      std::vector<ftype> v;
      v.resize(2 * n, 0);

      unpack_complex(in, v);

      gsl_fft_complex_wavetable *wave;
      gsl_fft_complex_workspace *work;

      work = gsl_fft_complex_workspace_alloc(n);
      wave = gsl_fft_complex_wavetable_alloc(n);

      gsl_fft_complex_inverse(v.data(), 1, n, wave, work);


      out.clear();

      pack_to_complex(v, out);

      out.resize(n, 0);

      gsl_fft_complex_wavetable_free(wave);
      gsl_fft_complex_workspace_free(work);
   }


   // Inverse of rfft
   // @in: input vector which must be the result of a rfft
   // @out: irfft of input, always real
   // Missing n: size of output
   // TODO fix this one!!
   static inline void irfft(complex_vector_t in, f_vector_t &out, uint n = 0)
   {
      assert(in.size() > 1);

      uint last = in.size() - 2;
      n = (n == 0) ? 2 * (in.size() - 1) : n;

      if (n == 2 * in.size() - 1) {
         last = in.size() - 1;
      } else if (n == 2 * (in.size() - 1)) {
         ;
      } else {
         std::cerr << "[mymath::ifft] Size not supported!\n"
                   << "[mymath::ifft] default size: " << n
                   << " will be used\n";
      }

      for (uint i = last; i > 0; --i) {
         in.push_back(std::conj(in[i]));
      }

      complex_vector_t temp;
      mymath::ifft(in, n, temp);

      mymath::complex_to_real(temp, out);

   }


// Same as python's numpy.fft.rfftfreq
// @ n: window length
// @ d (optional) : Sample spacing
// @return: A vector of length (n div 2) + 1 of the sample frequencies
   static inline f_vector_t rfftfreq(const uint n, const ftype d = 1.0)
   {
      f_vector_t v(n / 2 + 1);
      const ftype factor = 1.0 / (d * n);
      #pragma omp parallel for
      for (uint i = 0; i < v.size(); ++i) {
         v[i] = i * factor;
      }
      return std::move(v);
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
