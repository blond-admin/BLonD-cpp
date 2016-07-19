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

class API  PhaseLoop {
  public:
    virtual void track(){};
    void default_track();
    PhaseLoop(f_vector_t PL_gain, ftype window_coefficient, uint _delay,
              PhaseNoise* phaseNoise, LHCNoiseFB* LHCNoiseFB);
    void beam_phase();
    void phase_difference();
    void radial_steering_from_freq();
    void radial_difference();
    PhaseLoop(){};
    uint delay = 0;
    ftype alpha = 0;
    // f_vector_t gain;
    // f_vector_t gain2;
    ftype domega_RF = 0;
    ftype drho = 0;
    // f_vector_t domega_RF;
    ftype phi_beam = 0;
    ftype dphi = 0;
    ftype reference = 0;
    PhaseNoise* RFnoise;
    LHCNoiseFB* noiseFB;
    virtual ~PhaseLoop(){};
};

class API  LHC : public PhaseLoop {
  private:
  public:
    ftype gain2;
    f_vector_t gain;
    ftype lhc_y;
    f_vector_t lhc_a;
    f_vector_t lhc_t;

    ~LHC();
    void track();
    LHC(f_vector_t PL_gain, ftype SL_gain = 0, ftype window_coefficient = 0,
        PhaseNoise* phaseNoise = NULL, LHCNoiseFB* LHCNoiseFB = NULL,
        uint _delay = 0);
};

class API  PSB : public PhaseLoop {
  private:
  public:
    f_vector_t gain2;
    f_vector_t gain;
    f_vector_t dt;
    uint PL_counter;
    uint_vector_t on_time;
    f_vector_t coefficients;
    ftype dphi_av;
    ftype dphi_av_prev;
    ftype drho_prev;
    ftype t_accum;
    ftype domega_PL;
    ftype domega_RL;
    ~PSB();
    void track();
    PSB(f_vector_t PL_gain, f_vector_t RL_gain = f_vector_t(),
        ftype PL_period = 0, ftype RL_period = 0,
        f_vector_t coefficients = f_vector_t(), ftype window_coefficient = 0,
        PhaseNoise* phaseNoise = NULL, LHCNoiseFB* LHCNoiseFB = NULL,
        uint delay = 0);
    void precalculate_time();
};

class API  LHC_F : public PhaseLoop {
  private:
  public:
    ftype gain2;
    // f_vector_t domega_RF;
    ftype gain;
    LHC_F(ftype PL_gain, ftype window_coefficient = 0, ftype FL_gain = 0,
          PhaseNoise* phaseNoise = NULL, LHCNoiseFB* LHCNoiseFB = NULL,
          uint delay = 0);
    ~LHC_F();
    void track();
};

class API  SPS_RL : public PhaseLoop {
  private:
  public:
    ftype gain2;
    // f_vector_t domega_RF;
    ftype gain;

    SPS_RL(ftype PL_gain, ftype window_coefficient = 0, ftype RL_gain = 0,
           PhaseNoise* phaseNoise = NULL, LHCNoiseFB* LHCNoiseFB = NULL,
           uint delay = 0);
    ~SPS_RL();
    void track();
};

#endif /* PHASELOOP_H_ */
