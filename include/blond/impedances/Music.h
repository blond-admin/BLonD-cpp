/*
 * Music.h
 *
 *  Created on: Nov 25, 2016
 *      Author: kiliakis
 */

#pragma once

#include <blond/beams/Beams.h>
#include <blond/configuration.h>
#include <blond/impedances/Intensity.h>
namespace blond {
    class Music {
    private:

    public:
        Beams *Beam;
        double R_S;
        double omega_R;
        double Q;
        int n_macroparticles;
        long int n_particles;
        double alpha;
        double omega_bar;
        double constant;
        f_vector_t induced_voltage;
        double coeff1;
        double coeff2;
        double coeff3;
        double coeff4;

        Music(Beams *beam, const f_vector_t &resonator,
              int n_macroparticles, long int n_particles);
        ~Music() {};
        void track();
        void track_classic();
    };

} // blond