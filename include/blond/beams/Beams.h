/*
 * Beams.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_BEAMS_H_
#define BEAMS_BEAMS_H_

class Beams;

#include <blond/configuration.h>
#include <blond/utilities.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>

class API Beams {
public:
    f_vector_t dt;
    f_vector_t dE;
    // #NOTE id is 1 for active, 0 for inactive particles
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
    void losses_longitudinal_cut(const ftype dt_min, const ftype dt_max);
    void losses_energy_cut(const ftype dE_min, const ftype dE_max);
    void losses_separatrix(GeneralParameters *GP, RfParameters *RfP);

    // void losses_longitudinal_cut(const ftype* __restrict dt, const ftype dt_min,
    //                              const ftype dt_max, int* __restrict id);
    // void losses_energy_cut(const ftype* __restrict dE, const ftype dE_min,
    //                        const ftype dE_max, int* __restrict id);

    void statistics();

private:
    void statistics(const double *__restrict dE,
                    const double *__restrict dt,
                    const int *__restrict id,
                    const int size);
};

#endif /* BEAMS_BEAMS_H_ */
