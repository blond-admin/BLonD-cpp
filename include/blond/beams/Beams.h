/*
 * Beams.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_BEAMS_H_
#define BEAMS_BEAMS_H_

#include <blond/configuration.h>
#include <blond/utilities.h>

class API Beams {
  public:
    f_vector_t dt;
    f_vector_t dE;
    int_vector_t id;

    ftype mean_dt;
    ftype mean_dE;
    ftype sigma_dt;
    ftype sigma_dE;
    ftype ratio;
    ftype epsn_rms_l;
    uint n_macroparticles_lost;
    uint n_macroparticles;
    long long intensity;
    Beams(const uint _n_macroparticles, const long long _intensity);
    ~Beams();
    uint n_macroparticles_alive();
    void losses_longitudinal_cut(const ftype* __restrict dt, const ftype dt_min,
                                 const ftype dt_max, int* __restrict id);
    void losses_energy_cut(const ftype* __restrict dE, const ftype dE_min,
                           const ftype dE_max, int* __restrict id);

    void statistics();
};

#endif /* BEAMS_BEAMS_H_ */
