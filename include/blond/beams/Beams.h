/*
 * Beams.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_BEAMS_H_
#define BEAMS_BEAMS_H_

class Beams;

//#include <blond/input_parameters/GeneralParameters.h>
#include <blond/globals.h>
#include <blond/utilities.h>
#include <blond/configuration.h>
#include <blond/constants.h>

class Beams {
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
    int n_macroparticles_lost;
    int n_macroparticles;
    long intensity;

    //GeneralParameters *gp;
    Beams(const int _n_macroparticles, const long _intensity);

    ~Beams();

    int n_macroparticles_alive();

    void losses_longitudinal_cut(const ftype *__restrict__ dt,
                                 const ftype dt_min, const ftype dt_max, int *__restrict__ id);

    void losses_energy_cut(const ftype *__restrict__ dE, const ftype dE_min,
                           const ftype dE_max, int *__restrict__ id);

    void statistics();
};

#endif /* BEAMS_BEAMS_H_ */
