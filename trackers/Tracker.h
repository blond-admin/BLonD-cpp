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
   void set_periodicity();
   // Periodicity kick
   void kick(const int_vector_t &filter, const int index);
   inline void kick(const ftype *__restrict__ beam_dt,
                    ftype *__restrict__ beam_dE, const int n_rf,
                    const ftype *__restrict__ voltage,
                    const ftype *__restrict__ omega_RF,
                    const ftype *__restrict__ phi_RF, const int n_macroparticles,
                    const ftype acc_kick, const int_vector_t &filter);
   // Regular kick
   inline void kick(const int index);
   inline void kick(const ftype *__restrict__ beam_dt,
                    ftype *__restrict__ beam_dE, const int n_rf,
                    const ftype *__restrict__ voltage,
                    const ftype *__restrict__ omega_RF,
                    const ftype *__restrict__ phi_RF, const int n_macroparticles,
                    const ftype acc_kick);

   // Periodicity drift
   void drift(const int_vector_t &filter, const int index);
   inline void drift(ftype *__restrict__ beam_dt,
                     const ftype *__restrict__ beam_dE, const solver_type solver,
                     const ftype T0, const ftype length_ratio, const int alpha_order,
                     const ftype eta_zero, const ftype eta_one, const ftype eta_two,
                     const ftype beta, const ftype energy, const int n_macroparticles,
                     const int_vector_t &filter);
   // Regular drift
   inline void drift(const int index);
   inline void drift(ftype *__restrict__ beam_dt,
                     const ftype *__restrict__ beam_dE, const solver_type solver,
                     const ftype T0, const ftype length_ratio, const int alpha_order,
                     const ftype eta_zero, const ftype eta_one, const ftype eta_two,
                     const ftype beta, const ftype energy, const int n_macroparticles);

   //void track(const int start, const int end);

   void track();

   inline void horizontal_cut();
   RingAndRfSection(solver_type solver = simple, PhaseLoop *PL = NULL,
                    ftype *NoiseFB = NULL, bool periodicity = false, ftype dE_max = 0,
                    bool rf_kick_interp = false, ftype *Slices = NULL,
                    ftype *TotalInducedVoltage = NULL);
   ~RingAndRfSection();

   solver_type solver;
   PhaseLoop *PL;
   ftype *NoiseFB;
   bool periodicity;
   ftype dE_max;
   bool rf_kick_interp;
   ftype *Slices;
   ftype *TotalInducedVoltage;
   //int n_threads;
   ftype *acceleration_kick;

};

#endif /* TRACKERS_TRACKER_H_ */
