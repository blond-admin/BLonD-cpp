/*
 * PhaseLoop.h
 *
 *  Created on: Apr 7, 2016
 *      Author: kiliakis
 */

#ifndef PHASELOOP_H_
#define PHASELOOP_H_


#include "globals.h"
#include "utilities.h"
#include "configuration.h"
#include "constants.h"
#include "PhaseNoise.h"
#include "LHCNoiseFB.h"


class PhaseLoop {
public:
   virtual void track() {};
   void default_track();
   PhaseLoop(f_vector_t PL_gain,
             ftype window_coefficient,
             uint _delay,
             PhaseNoise *phaseNoise,
             LHCNoiseFB *LHCNoiseFB);
   void beam_phase();
   void phase_difference();
   void radial_steering_from_freq();
   PhaseLoop() {} ;
   uint delay = 0;
   ftype alpha = 0;
   f_vector_t gain;
   ftype drho = 0;
   ftype domega_RF = 0;
   ftype phi_beam = 0;
   ftype dphi = 0;
   ftype reference = 0;
   PhaseNoise *RFnoise;
   LHCNoiseFB *noiseFB;
   virtual ~PhaseLoop() {};
};

class LHC: public PhaseLoop {
private:
   //ftype gain;
   //ftype domega_RF;
public:
   ftype gain2;
   ftype lhc_y;
   f_vector_t lhc_a;
   f_vector_t lhc_t;

   ~LHC();
   void track();
   LHC(f_vector_t PL_gain,
       ftype SL_gain = 0,
       ftype window_coefficient = 0,
       PhaseNoise *phaseNoise = NULL,
       LHCNoiseFB *LHCNoiseFB = NULL,
       uint _delay = 0);
};

class PSB: public PhaseLoop {
private:
   f_vector_t gain2;
   //f_vector_t gain;
   uint_vector_t dt;
   ftype average_dE;
   uint PL_counter;
   uint_vector_t on_time;
   f_vector_t coefficients;
   ftype dphi_av;
   ftype dphi_av_prev;
   ftype drho_prev;
   ftype t_accum;
   ftype domega_PL;
   ftype domega_RL;
   //ftype domega_RF;
public:
   ~PSB();
   void track();
   PSB(f_vector_t PL_gain,
       f_vector_t RL_gain = f_vector_t(),
       ftype PL_period = 0,
       ftype RL_period = 0,
       f_vector_t coefficients = f_vector_t(),
       ftype window_coefficient = 0,
       PhaseNoise *phaseNoise = NULL,
       LHCNoiseFB *LHCNoiseFB = NULL,
       uint delay = 0);
   void radial_difference();
   void precalculate_time();
};

#endif /* PHASELOOP_H_ */
