/*
 * PhaseLoop.h
 *
 *  Created on: Apr 7, 2016
 *      Author: kiliakis
 */

#ifndef PHASELOOP_H_
#define PHASELOOP_H_

#include <blond/configuration.h>
#include <blond/globals.h>
#include <blond/llrf/LHCNoiseFB.h>
#include <blond/llrf/PhaseNoise.h>
#include <blond/utilities.h>
namespace blond {
    class PhaseLoop {
    public:
        virtual void track() {};
        void default_track();
        PhaseLoop(f_vector_t PL_gain, double window_coefficient, uint _delay,
                  PhaseNoise *phaseNoise, LHCNoiseFB *LHCNoiseFB);
        void beam_phase();
        void phase_difference();
        void radial_steering_from_freq();
        void radial_difference();
        PhaseLoop() {};
        uint delay = 0;
        double alpha = 0;
        // f_vector_t gain;
        // f_vector_t gain2;
        double domega_rf = 0;
        double drho = 0;
        // f_vector_t domega_RF;
        double phi_beam = 0;
        double dphi = 0;
        double reference = 0;
        PhaseNoise *RFnoise;
        LHCNoiseFB *noiseFB;
        virtual ~PhaseLoop() {};
    };

    class LHC : public PhaseLoop {
    private:
    public:
        double gain2;
        f_vector_t gain;
        double lhc_y;
        f_vector_t lhc_a;
        f_vector_t lhc_t;

        ~LHC();
        void track();
        LHC(f_vector_t PL_gain, double SL_gain = 0, double window_coefficient = 0,
            PhaseNoise *phaseNoise = NULL, LHCNoiseFB *LHCNoiseFB = NULL,
            uint _delay = 0);
    };

    class PSB : public PhaseLoop {
    private:
    public:
        f_vector_t gain2;
        f_vector_t gain;
        f_vector_t dt;
        uint PL_counter;
        uint_vector_t on_time;
        f_vector_t coefficients;
        double dphi_av;
        double dphi_av_prev;
        double drho_prev;
        double t_accum;
        double domega_PL;
        double domega_RL;
        ~PSB();
        void track();
        PSB(f_vector_t PL_gain, f_vector_t RL_gain = f_vector_t(),
            double PL_period = 0, double RL_period = 0,
            f_vector_t coefficients = f_vector_t(), double window_coefficient = 0,
            PhaseNoise *phaseNoise = NULL, LHCNoiseFB *LHCNoiseFB = NULL,
            uint delay = 0);
        void precalculate_time();
    };

    class LHC_F : public PhaseLoop {
    private:
    public:
        double gain2;
        // f_vector_t domega_RF;
        double gain;
        LHC_F(double PL_gain, double window_coefficient = 0, double FL_gain = 0,
              PhaseNoise *phaseNoise = NULL, LHCNoiseFB *LHCNoiseFB = NULL,
              uint delay = 0);
        ~LHC_F();
        void track();
    };

    class SPS_RL : public PhaseLoop {
    private:
    public:
        double gain2;
        // f_vector_t domega_RF;
        double gain;

        SPS_RL(double PL_gain, double window_coefficient = 0, double RL_gain = 0,
               PhaseNoise *phaseNoise = NULL, LHCNoiseFB *LHCNoiseFB = NULL,
               uint delay = 0);
        ~SPS_RL();
        void track();
    };
} // blond
#endif /* PHASELOOP_H_ */
