/*
 * math_functions.h
 *
 *  Created on: Mar 21, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDE_FFT_H_
#define INCLUDE_FFT_H_

#include <algorithm>
#include <blond/configuration.h>
#include <blond/utilities.h>
#include <cassert>
#include <cmath>
#include <fftw3.h>
#include <functional>
#include <blond/openmp.h>

namespace fft {

    // FFTW_PATIENT: run a lot of ffts to discover the best plan.
    // Will not use the multithreaded version unless the fft size
    // is reasonable enough

    // FFTW_MEASURE : run some ffts to find the best plan.
    // (first run should take some more seconds)

    // FFTW_ESTIMATE : dont run any ffts, just make an estimation.
    // Usually leads to suboptimal solutions

    // FFTW_DESTROY_INPUT : use the original input to store arbitaty data.
    // May yield better performance but the input is not usable any more.
    // Can be combined with all the above
    const uint FFTW_FLAGS = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;
    // const uint FFTW_FLAGS = FFTW_ESTIMATE;// | FFTW_DESTROY_INPUT;
    // const uint ELEMS_PER_THREAD_FFT = 10000;
    // const uint ELEMS_PER_THREAD_RFFT = 15000;

    enum fft_type_t { FFT, IFFT, RFFT, IRFFT };

    struct API fft_plan_t {
        fftw_plan p; // fftw_plan
        uint n;      // size of the fft
        fft_type_t type;
        void *in;
        void *out;
    };

    static std::vector<fft::fft_plan_t> planV;



    static inline void real_to_complex(const std::vector<double> &in,
                                       std::vector<complex_t> &out)
    {
        assert(out.empty());
        out.reserve(in.size());
        for (const auto &real : in)
            out.push_back(complex_t(real, 0));
    }

    static inline void pack_to_complex(const std::vector<double> &in,
                                       std::vector<complex_t> &out)
    {
        assert(out.empty());
        assert(in.size() % 2 == 0);
        out.reserve(in.size() / 2);
        for (unsigned int i = 0; i < in.size(); i += 2)
            out.push_back(complex_t(in[i], in[i + 1]));
    }

    static inline void complex_to_real(const std::vector<complex_t> &in,
                                       std::vector<double> &out)
    {
        assert(out.empty());
        out.reserve(in.size());
        for (const auto &z : in)
            out.push_back(z.real());
    }

    static inline void unpack_complex(const std::vector<complex_t> &in,
                                      std::vector<double> &out)
    {
        assert(out.empty());
        out.reserve(2 * in.size());
        for (const auto &z : in) {
            out.push_back(z.real());
            out.push_back(z.imag());
        }
    }


    //#ifdef USE_FFTW
    static inline fftw_plan init_fft(const int n, complex_t *in, complex_t *out,
                                     const int sign = FFTW_FORWARD,
                                     const unsigned flag = FFTW_ESTIMATE,
                                     const int threads = 1)
    {
#ifdef USE_FFTW_OMP
        if (threads > 1) {
            fftw_init_threads();
            fftw_plan_with_nthreads(threads);
            // std::min(
            //     threads,
            //     (int)((n + ELEMS_PER_THREAD_FFT - 1) / ELEMS_PER_THREAD_FFT)));
        }
#endif
        fftw_complex *a, *b;
        a = reinterpret_cast<fftw_complex *>(in);
        b = reinterpret_cast<fftw_complex *>(out);
        return fftw_plan_dft_1d(n, a, b, sign, flag);
    }

    static inline fftw_plan init_rfft(const int n, double *in, complex_t *out,
                                      const unsigned flag = FFTW_ESTIMATE,
                                      const int threads = 1)

    {
#ifdef USE_FFTW_OMP
        if (threads > 1) {
            fftw_init_threads();
            fftw_plan_with_nthreads(threads);
            // std::min(threads, (int)((n + ELEMS_PER_THREAD_RFFT - 1) /
            //                         ELEMS_PER_THREAD_RFFT)));
        }
#endif
        fftw_complex *b;
        b = reinterpret_cast<fftw_complex *>(out);
        return fftw_plan_dft_r2c_1d(n, in, b, flag);
    }

    static inline fftw_plan init_irfft(const int n, complex_t *in, double *out,
                                       const unsigned flag = FFTW_ESTIMATE,
                                       const int threads = 1)
    {
#ifdef USE_FFTW_OMP
        if (threads > 1) {
            fftw_init_threads();
            fftw_plan_with_nthreads(threads);
            // std::min(threads, (int)((n + ELEMS_PER_THREAD_FFT - 1) /
            //                         ELEMS_PER_THREAD_FFT)));
        }
#endif
        fftw_complex *b;
        b = reinterpret_cast<fftw_complex *>(in);
        return fftw_plan_dft_c2r_1d(n, b, out, flag);
    }

    static inline void run_fft(const fftw_plan &p) { fftw_execute(p); }

    static inline void destroy_fft(fftw_plan &p) { fftw_destroy_plan(p); }

    static inline void destroy_plans()
    {
        // std::cout<<"planV size = " << planV.size()<<"\n";
        for (auto &i : planV) {
            fftw_destroy_plan(i.p);
            fftw_free(i.in);
            fftw_free(i.out);
        }
        planV.clear();
    }

    //#endif

    static inline fft_plan_t find_plan(uint n, fft_type_t type, uint threads,
                                       std::vector<fft_plan_t> &v)
    {
        const uint flag = FFTW_FLAGS;
        auto it =
        std::find_if(v.begin(), v.end(), [n, type](const fft_plan_t &s) {
            return ((s.n == n) && (s.type == type));
        });

        if (it == v.end()) {
            // std::cout << "I have to create a new plan :(\n";
            fft_plan_t plan;
            plan.n = n;
            plan.type = type;

            if (type == FFT) {
                fftw_complex *in =
                    (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
                fftw_complex *out =
                    (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);

                auto p = init_fft(n, reinterpret_cast<complex_t *>(in),
                                  reinterpret_cast<complex_t *>(out),
                                  FFTW_FORWARD, flag, threads);
                plan.p = p;
                plan.in = in;
                plan.out = out;
            } else if (type == IFFT) {

                fftw_complex *in =
                    (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
                fftw_complex *out =
                    (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
                auto p = init_fft(n, reinterpret_cast<complex_t *>(in),
                                  reinterpret_cast<complex_t *>(out),
                                  FFTW_BACKWARD, flag, threads);
                plan.p = p;
                plan.in = in;
                plan.out = out;

            } else if (type == RFFT) {
                double *in = (double *)fftw_malloc(sizeof(double) * n);
                fftw_complex *out = (fftw_complex *)fftw_malloc(
                                        sizeof(fftw_complex) * (n / 2 + 1));
                auto p = init_rfft(n, in, reinterpret_cast<complex_t *>(out),
                                   flag, threads);
                // std::cout << "threads : " << threads << "\n";
                plan.p = p;
                plan.in = in;
                plan.out = out;

            } else if (type == IRFFT) {
                double *out = (double *)fftw_malloc(sizeof(double) * n);
                fftw_complex *in = (fftw_complex *)fftw_malloc(
                                       sizeof(fftw_complex) * (n / 2 + 1));
                auto p = init_irfft(n, reinterpret_cast<complex_t *>(in), out,
                                    flag, threads);
                plan.p = p;
                plan.in = in;
                plan.out = out;

            } else {
                std::cerr << "[fft::find_plan]: ERROR "
                          << "Wrong fft type!\n";
                exit(-1);
            }

            v.push_back(plan);
            // planV.push_back(plan);
            return plan;
        } else {
            return *it;
        }
    }

    // Parameters are like python's numpy.fft.rfft
    // @in:  input data
    // @n:   number of points to use. If n < in.size() then the input is cropped
    //       if n > in.size() then input is padded with zeros
    // @out: the transformed array
    static inline void rfft(f_vector_t &in, complex_vector_t &out, uint n = 0,
                            const uint threads = 1)
    {
        if (n == 0)
            n = in.size();
        else
            in.resize(n);

        out.resize(n / 2 + 1);

        auto plan = fft::find_plan(n, RFFT, threads, planV);
        auto from = (double *)plan.in;
        auto to = (complex_t *)plan.out;

        std::copy(in.begin(), in.end(), from);

        run_fft(plan.p);

        std::copy(&to[0], &to[n / 2 + 1], out.begin());
    }

    // Parameters are like python's numpy.fft.fft
    // @in:  input data
    // @n:   number of points to use. If n < in.size() then the input is cropped
    //       if n > in.size() then input is padded with zeros
    // @out: the transformed array
    static inline void fft(complex_vector_t &in, complex_vector_t &out,
                           uint n = 0, const uint threads = 1)
    {
        if (n == 0)
            n = in.size();

        out.resize(n);

        auto plan = fft::find_plan(n, FFT, threads, planV);
        auto from = (complex_t *)plan.in;
        auto to = (complex_t *)plan.out;

        std::copy(in.begin(), in.end(), from);

        run_fft(plan.p);

        std::copy(&to[0], &to[n], out.begin());
    }

    // Parameters are like python's numpy.fft.ifft
    // @in:  input data
    // @n:   number of points to use. If n < in.size() then the input is cropped
    //       if n > in.size() then input is padded with zeros
    // @out: the inverse Fourier transform of input data
    static inline void ifft(complex_vector_t &in, complex_vector_t &out,
                            uint n = 0, const uint threads = 1)
    {
        if (n == 0)
            n = in.size();

        out.resize(n);

        auto plan = fft::find_plan(n, IFFT, threads, planV);
        auto from = (complex_t *)plan.in;
        auto to = (complex_t *)plan.out;
        std::copy(in.begin(), in.end(), from);

        run_fft(plan.p);

        std::transform(&to[0], &to[n], out.begin(),
                       std::bind2nd(std::divides<complex_t>(), n));
    }

    // Inverse of rfft
    // @in: input vector which must be the result of a rfft
    // @out: irfft of input, always real
    // Missing n: size of output
    static inline void irfft(complex_vector_t in, f_vector_t &out, uint n = 0,
                             const uint threads = 1)
    {
        n = (n == 0) ? 2 * (in.size() - 1) : n;
        // std::cout << "out size will be " << n << "\n";
        out.resize(n);

        auto plan = fft::find_plan(n, IRFFT, threads, planV);
        auto from = (complex_t *)plan.in;
        auto to = (double *)plan.out;

        std::copy(in.begin(), in.end(), from);

        run_fft(plan.p);

        std::transform(&to[0], &to[n], out.begin(),
                       std::bind2nd(std::divides<double>(), n));
    }

    // Same as python's numpy.fft.rfftfreq
    // @ n: window length
    // @ d (optional) : Sample spacing
    // @return: A vector of length (n div 2) + 1 of the sample frequencies
    static inline f_vector_t rfftfreq(const uint n, const double d = 1.0)
    {
        f_vector_t v(n / 2 + 1);
        const double factor = 1.0 / (d * n);
        #pragma omp parallel for
        for (int i = 0; i < (int)v.size(); ++i) {
            v[i] = i * factor;
        }
        return v;
    }


    static inline void convolution_with_ffts(f_vector_t signal,
            f_vector_t kernel,
            f_vector_t &res)
    {
        complex_vector_t v1; //(signal.size());
        complex_vector_t v2; //(kernel.size());
        const long unsigned size = signal.size() + kernel.size() - 1;
        res.resize(size);

        fft::rfft(signal, v1, size, omp_get_max_threads());
        fft::rfft(kernel, v2, size, omp_get_max_threads());

        std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(),
                       std::multiplies<complex_t>());

        fft::irfft(v1, res, size, omp_get_max_threads());
    }
}

#endif /* INCLUDE_FFT_H_ */
