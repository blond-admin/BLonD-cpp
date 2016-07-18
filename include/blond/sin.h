/*
 * sin.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#include <blond/sincos.h>

#ifndef TRACKERS_SIN_H_
#define TRACKERS_SIN_H_

namespace vdt {

    // Sin double precision
    // --------------------------------------------------------

    /// Double precision sine: just call sincos.
    inline double fast_sin(double x) {
        double s, c;
        fast_sincos(x, s, c);
        return s;
    }

    //------------------------------------------------------------------------------

    inline float fast_sinf(float x) {
        float s, c;
        fast_sincosf(x, s, c);
        return s;
    }

    //------------------------------------------------------------------------------
    void sinv(const uint32_t size, double const* __restrict iarray,
              double* __restrict oarray);

    void fast_sinv(const uint32_t size, double const* __restrict iarray,
                   double* __restrict oarray);

    void sinfv(const uint32_t size, float const* __restrict iarray,
               float* __restrict oarray);

    void fast_sinfv(const uint32_t size, float const* __restrict iarray,
                    float* __restrict oarray);

} // vdt namespace

#endif /* TRACKERS_SIN_H_ */
