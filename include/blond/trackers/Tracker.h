/*
 * Tracker.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#ifndef TRACKERS_TRACKER_H_
#define TRACKERS_TRACKER_H_

#include "utilities.h"
#include "../input_parameters/GeneralParameters.h"
#include "../input_parameters/RfParameters.h"
#include "../beams/Beams.h"
#include "../llrf/PhaseLoop.h"
#include "../llrf/LHCNoiseFB.h"
#include "../beams/Slices.h"
#include "../impedances/InducedVoltage.h"

enum solver_type {
   simple, full
};


class RingAndRfSection {

private:
public:
   ftype elapsed_time;
   int_vector_t indices_right_outside;
   int_vector_t indices_inside_frame;
   int_vector_t indices_left_outside;

   f_vector_t acceleration_kick;
   solver_type solver;
   ftype dE_max;
   bool rf_kick_interp;
   bool periodicity;

   LHCNoiseFB *noiseFB;
   PhaseLoop *PL;
   Slices *slices;
   TotalInducedVoltage *totalInducedVoltage;

   //helping arrays for kick
   // ftype *vol, *omeg, *phi; 

   void set_periodicity();
   // Periodicity kick
   void kick(const int_vector_t &filter,
             const uint index);
   inline void kick(const ftype *__restrict__ beam_dt,
                    ftype *__restrict__ beam_dE,
                    const int n_rf,
                    const ftype *__restrict__ voltage,
                    const ftype *__restrict__ omega_RF,
                    const ftype *__restrict__ phi_RF,
                    const int n_macroparticles,
                    const ftype acc_kick,
                    const int_vector_t &filter);
   // Regular kick
   inline void kick(const uint index);
   inline void kick(const ftype *__restrict__ beam_dt,
                    ftype *__restrict__ beam_dE,
                    const int n_rf,
                    const ftype *__restrict__ voltage,
                    const ftype *__restrict__ omega_RF,
                    const ftype *__restrict__ phi_RF,
                    const int n_macroparticles,
                    const ftype acc_kick);

   // Periodicity drift
   void drift(const int_vector_t &filter,
              const uint index);
   inline void drift(ftype *__restrict__ beam_dt,
                     const ftype *__restrict__ beam_dE,
                     const solver_type solver,
                     const ftype T0,
                     const ftype length_ratio,
                     const uint alpha_order,
                     const ftype eta_zero,
                     const ftype eta_one,
                     const ftype eta_two,
                     const ftype beta,
                     const ftype energy,
                     const int n_macroparticles,
                     const int_vector_t &filter);
   // Regular drift
   inline void drift(const uint index);
   inline void drift(ftype *__restrict__ beam_dt,
                     const ftype *__restrict__ beam_dE,
                     const solver_type solver,
                     const ftype T0,
                     const ftype length_ratio,
                     const uint alpha_order,
                     const ftype eta_zero,
                     const ftype eta_one,
                     const ftype eta_two,
                     const ftype beta,
                     const ftype energy,
                     const int n_macroparticles);


   void track();

   inline void horizontal_cut();
   RingAndRfSection(solver_type solver = simple,
                    PhaseLoop *PL = NULL,
                    LHCNoiseFB *NoiseFB = NULL,
                    bool periodicity = false,
                    ftype dE_max = 0,
                    bool rf_kick_interp = false,
                    Slices *Slices = NULL,
                    TotalInducedVoltage *TotalInducedVoltage = NULL);
   ~RingAndRfSection();

};

#endif /* TRACKERS_TRACKER_H_ */
