/*
 * Beams.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_BEAMS_H_
#define BEAMS_BEAMS_H_

class Beams;

#include <blond/globals.h>
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
    double mass;
    double charge;
    double beta;
    double gamma;
    double energy;
    double momentum;
    double mean_dt;
    double mean_dE;
    double sigma_dt;
    double sigma_dE;
    double ratio;
    double epsn_rms_l;
    int n_macroparticles_lost;
    int n_macroparticles;
    long long intensity;
    Beams(GeneralParameters *GP,
          const int _n_macroparticles,
          const long long _intensity);
    ~Beams();
    uint n_macroparticles_alive();
    void losses_longitudinal_cut(const double dt_min, const double dt_max);
    void losses_energy_cut(const double dE_min, const double dE_max);
    void losses_separatrix(GeneralParameters *GP, RfParameters *RfP);

    // void losses_longitudinal_cut(const double* __restrict dt, const double dt_min,
    //                              const double dt_max, int* __restrict id);
    // void losses_energy_cut(const double* __restrict dE, const double dE_min,
    //                        const double dE_max, int* __restrict id);

    void statistics();

private:
    void statistics(const double *__restrict dE,
                    const double *__restrict dt,
                    const int *__restrict id,
                    const int size);
};

#endif /* BEAMS_BEAMS_H_ */
