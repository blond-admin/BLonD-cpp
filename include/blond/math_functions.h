/*
 * math_functions.h
 *
 *  Created on: Mar 21, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_MATH_FUNCTIONS_H_
#define INCLUDES_MATH_FUNCTIONS_H_

#include <algorithm>
#include <blond/configuration.h>
#include <blond/exp.h>
#include <blond/sin.h>
#include <blond/utilities.h>
#include <cmath>
#include <random>
#include <limits>
#include <cassert>
#include <blond/openmp.h>

namespace mymath {

    class API Ham {
    private:
        std::vector<uint> _H, _hp, _hv, _x;

    public:
        bool operator!=(const Ham &other) const { return true; }

        Ham begin() const { return *this; }

        Ham end() const { return *this; }

        uint operator*() const { return _x.back(); }

        Ham(const std::vector<uint> &pfs)
            : _H(pfs), _hp(pfs.size(), 0), _hv(pfs), _x(1, 1) {}

        const Ham &operator++()
        {
            for (uint i = 0; i < _H.size(); i++)
                for (; _hv[i] <= _x.back(); _hv[i] = _x[++_hp[i]] * _H[i])
                    ;
            _x.push_back(_hv[0]);
            for (uint i = 1; i < _H.size(); i++)
                if (_hv[i] < _x.back())
                    _x.back() = _hv[i];
            return *this;
        }
    };

    static inline uint next_regular(uint target)
    {

        for (auto i : Ham({2, 3, 5})) {
            if (i > target)
                return i;
        }
        return 0;
    }

    // Wrapper function for vdt::fucntions
    static inline double fast_sin(double x) { return vdt::fast_sin(x); }

    static inline double fast_cos(double x) { return vdt::fast_sin(x + M_PI_2); }

    static inline double fast_exp(double x) { return vdt::fast_exp(x); }

    // linear convolution function
    static inline void convolution(const double *__restrict signal,
                                   const int SignalLen,
                                   const double *__restrict kernel,
                                   const int KernelLen, double *__restrict res)
    {
        const int size = KernelLen + SignalLen - 1;

        #pragma omp parallel for
        for (auto n = 0; n < size; ++n) {
            res[n] = 0;
            const int kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
            const int kmax = (n < SignalLen - 1) ? n : SignalLen - 1;
            for (int k = kmin; k <= kmax; k++) {
                res[n] += signal[k] * kernel[n - k];
                //--j;
            }
        }
    }


    // Parameters are like python's np.interp
    // @x: x-coordinates of the interpolated values
    // @xp: The x-coords of the data points
    // @fp: the y-coords of the data points
    // @y: the interpolated values, same shape as x
    // @left: value to return for x < xp[0]
    // @right: value to return for x > xp[last]
    static inline void interp(const std::vector<double> &x,
                              const std::vector<double> &xp,
                              const std::vector<double> &yp,
                              std::vector<double> &y,
                              const double left,
                              const double right)
    {
        assert(xp.size() == yp.size());
        assert(std::is_sorted(xp.begin(), xp.end()));
        y.resize(x.size());

        const int M = xp.size();
        const int N = x.size();
        if (N == 0) return;

        // const auto max = xp.back();
        // const auto min = xp.front();

        // #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            int pos = std::lower_bound(xp.begin(), xp.end(), x[i]) - xp.begin();
            if (pos == M)
                y[i] = right;
            else if (xp[pos] == x[i])
                y[i] = yp[pos];
            else if (pos == 0)
                y[i] = left;
            else {
                y[i] = yp[pos - 1] +
                       (yp[pos] - yp[pos - 1]) * (x[i] - xp[pos - 1]) /
                       (xp[pos] - xp[pos - 1]);
            }
        }
    }

    static inline void interp(const std::vector<double> &x,
                              const std::vector<double> &xp,
                              const std::vector<double> &yp,
                              std::vector<double> &y)
    {
        auto left = yp.front();
        auto right = yp.back();
        interp(x, xp, yp, y, left, right);
    }

    static inline f_vector_t interp(const f_vector_t &x,
                                    const f_vector_t &xp,
                                    const f_vector_t &yp)
    {
        f_vector_t res;
        interp(x, xp, yp, res);
        return res;
    }

    static inline f_vector_t interp(const f_vector_t &x,
                                    const f_vector_t &xp,
                                    const f_vector_t &yp,
                                    const double left,
                                    const double right)
    {
        f_vector_t res;
        interp(x, xp, yp, res, left, right);
        return res;
    }

    // Function to implement integration of f(x) over the interval
    // [a,b] using the trapezoid rule with nsub subdivisions.
    template <typename T>
    static inline f_vector_t cum_trapezoid(const T *f,
                                           const T deltaX,
                                           const int nsub)
    {
        // initialize the partial sum to be f(a)+f(b) and
        // deltaX to be the step size using nsub subdivisions
        // double *psum = new double[nsub];
        f_vector_t psum(nsub - 1);
        // psum[0] = 0;
        const auto half_dx = deltaX / 2.0;
        psum[0] = (f[1] + f[0]) * half_dx;
        // increment the partial sum
        // #pragma omp parallel for
        for (int i = 1; i < nsub - 1; ++i)
            psum[i] = psum[i - 1] + (f[i + 1] + f[i]) * half_dx;

        return psum;
    }

    // Function to implement integration of f(x) over the interval
    // [a,b] using the trapezoid rule with nsub subdivisions.
    template <typename T>
    static inline f_vector_t cum_trapezoid(const T *f,
                                           const T deltaX,
                                           const T initial,
                                           const int nsub)
    {
        // initialize the partial sum to be f(a)+f(b) and
        // deltaX to be the step size using nsub subdivisions
        // double *psum = new double[nsub];
        f_vector_t psum(nsub);
        const auto half_dx = deltaX / 2.0;

        psum[0] = initial;
        psum[1] = (f[1] + f[0]) * half_dx;
        // increment the partial sum
        for (int i = 2; i < nsub; ++i)
            psum[i] = psum[i - 1] + (f[i] + f[i - 1]) * half_dx;

        return psum;
    }

    template<typename T>
    static inline f_vector_t cum_trapezoid(const std::vector<T> &f,
                                           const T deltaX)
    {
        return cum_trapezoid<T>(f.data(), deltaX, f.size());
    }

    template<typename T>
    static inline f_vector_t cum_trapezoid(const std::vector<T> &f,
                                           const T deltaX,
                                           const T initial)
    {
        return cum_trapezoid<T>(f.data(), deltaX, initial, f.size());
    }


    template <typename T>
    static inline double trapezoid(const T *f, const double *deltaX, const int nsub)
    {
        // initialize the partial sum to be f(a)+f(b) and
        // deltaX to be the step size using nsub subdivisions

        double psum = 0.0;
// increment the partial sum
        #pragma omp parallel for reduction(+ : psum)
        for (int i = 1; i < nsub; ++i) {
            psum += (f[i] + f[i - 1]) * (deltaX[i] - deltaX[i - 1]);
        }

        // return approximation
        return psum / 2;
    }

    template <typename T>
    static inline double trapezoid(const std::vector<T> &f,
                                   const std::vector<double> &deltaX)
    {
        return trapezoid<T>(f.data(), deltaX.data(), f.size());
    }

    template <typename T>
    static inline double trapezoid(const T *f, const double deltaX, const int nsub)
    {
        // initialize the partial sum to be f(a)+f(b) and
        // deltaX to be the step size using nsub subdivisions
        double psum = f[0] + f[nsub - 1]; // f(a)+f(b);
// double deltaX = (b-a)/nsub;

// increment the partial sum
        #pragma omp parallel for reduction(+ : psum)
        for (int i = 1; i < nsub - 1; ++i) {
            psum += 2 * f[i];
        }

        // multiply the sum by the constant deltaX/2.0
        psum = (deltaX / 2) * psum;

        // return approximation
        return psum;
    }

    template <typename T>
    static inline double trapezoid(const std::vector<T> &f, const double deltaX)
    {
        return trapezoid<T>(f.data(), deltaX, f.size());
    }

    template <typename T>
    static inline uint min(T *a, uint size, uint step = 1)
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

    template <typename T>
    static inline uint max(T *a, uint size, uint step = 1)
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

    static inline void linspace(double *a, const double start, const double end,
                                const int n, const int keep_from = 0)
    {
        const double step = (end - start) / (n - 1);
        // double value = start;
        //#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            if (i >= keep_from)
                a[i - keep_from] = start + i * step;
            // value += step;
        }
    }

    static inline f_vector_t linspace(const double start, const double end,
                                      const int n)
    {
        f_vector_t res(n);
        linspace(res.data(), start, end, n);
        return res;
    }

    template <typename T>
    static inline std::vector<T> arange(T start, T stop, T step = 1)
    {
        std::vector<T> values;
        for (T value = start; value < stop; value += step)
            values.push_back(value);
        return values;
    }

    template <typename T>
    static inline double sum(const T data[], const int n)
    {
        double m = 0.0;
        #pragma omp parallel for reduction(+ : m)
        for (int i = 0; i < n; ++i) {
            m += data[i];
        }
        return m;
    }

    template <typename T>
    static inline double sum(const std::vector<T> &data)
    {
        return sum(data.data(), data.size());
    }

    template <typename T>
    static inline double mean(const T data[], const int n)
    {
        return sum(data, n) / n;
    }

    template <typename T>
    static inline double mean(const std::vector<T> data)
    {
        return mean(data.data(), data.size());
    }

    template <typename T>
    static inline double standard_deviation(const T data[], const int n,
                                            const double mean)
    {
        double sum_deviation = 0.0;
        #pragma omp parallel for reduction(+ : sum_deviation)
        for (int i = 0; i < n; ++i)
            sum_deviation += (data[i] - mean) * (data[i] - mean);
        return std::sqrt(sum_deviation / n);
    }

    template <typename T>
    static inline double standard_deviation(const std::vector<T> data,
                                            const double mean)
    {
        return standard_deviation(data.data(), data.size(), mean);
    }

    template <typename T>
    static inline double standard_deviation(const T data[], const int n)
    {
        const double mean = mymath::mean(data, n);
        double sum_deviation = 0.0;
        #pragma omp parallel for reduction(+ : sum_deviation)
        for (int i = 0; i < n; ++i)
            sum_deviation += (data[i] - mean) * (data[i] - mean);
        return std::sqrt(sum_deviation / n);
    }

    template <typename T>
    static inline double standard_deviation(const std::vector<T> data)
    {
        return standard_deviation(data.data(), data.size());
    }

    template <typename T> int sign(const T val)
    {
        return (T(0) < val) - (val < T(0));
    }


    template <typename T>
    static inline void meshgrid(const std::vector<T> &in1,
                                const std::vector<T> &in2,
                                std::vector< std::vector<T> > &out1,
                                std::vector< std::vector<T> > &out2)
    {
        out1.clear();
        out1.resize(in2.size(), in1);

        out2.clear();
        const int cols = in1.size();
        FOR(in2, e) out2.push_back(std::vector<T>(cols, *e));

    }

    template <typename T>
    static inline std::vector<T> flatten(const std::vector<std::vector<T> > &a)
    {
        std::vector<T> result;
        FOR(a, row) {
            result.insert(result.end(), (*row).begin(), (*row).end());
        }
        return result;
    }

    template <typename T>
    static inline std::vector<T> random_choice(const std::vector<T> &elems,
            const int size,
            std::vector<double> weights = std::vector<double>())
    {
        if (weights.empty()) weights.resize(elems.size(), 1);
        assert(weights.size() == elems.size());
        // calculate cumulative sum
        for (uint i = 1; i < weights.size(); i++) {
            weights[i] += weights[i - 1];
        }
        std::default_random_engine gen(0);
        std::uniform_real_distribution<double> d(0.0, weights.back());
        std::vector<T> result(size);
        for (int i = 0; i < size; i++) {
            auto w = d(gen);
            auto it = std::lower_bound(ALL(weights), w);
            result[i] = elems[it - weights.begin()];
        }
        return result;
    }

    static inline f_vector_t real_part(const std::vector<complex_t> &in)
    {
        f_vector_t res; res.reserve(in.size());
        for (const auto &z : in) res.push_back(z.real());
        return res;
    }

    static inline f_vector_t imag_part(const std::vector<complex_t> &in)
    {
        f_vector_t res; res.reserve(in.size());
        for (const auto &z : in) res.push_back(z.imag());
        return res;
    }
}
#endif /* INCLUDES_MATH_FUNCTIONS_H_ */
