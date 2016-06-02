/*
 * Tracker.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include "Tracker.h"

RingAndRfSection::~RingAndRfSection()
{
   util::delete_array(this->TotalInducedVoltage);
   util::delete_array(this->acceleration_kick);
   util::delete_array(this->Slices);

}

// Two versions of kick, drift one with periodicity and another without periodiciy
// First go the versions without periodicity

// Kick without periodicity
inline void RingAndRfSection::kick(const ftype *__restrict__ beam_dt,
                                   ftype *__restrict__ beam_dE, const int n_rf,
                                   const ftype *__restrict__ voltage, const ftype *__restrict__ omega_RF,
                                   const ftype *__restrict__ phi_RF, const int n_macroparticles,
                                   const ftype acc_kick)
{

   //beam_dE[0] += 1;
// KICK



   int k = 0;
   for (int j = 0; j < n_rf; j++) {
      //const ftype orf = omega_RF[k];
      //const ftype prf = phi_RF[k];
      //const ftype v = voltage[k];
      #pragma omp parallel for
      for (int i = 0; i < n_macroparticles; i++) {
         //beam_dE[i] += v * vdt::fast_sin(orf * beam_dt[i] + prf);
         beam_dE[i] += voltage[k] *
                       vdt::fast_sin(omega_RF[k] * beam_dt[i] + phi_RF[k]);
      }
      k += GP->n_turns;
   }

// SYNCHRONOUS ENERGY CHANGE
   #pragma omp parallel for
   for (int i = 0; i < n_macroparticles; i++)
      beam_dE[i] += acc_kick;

}

// kick with periodicity

inline void RingAndRfSection::kick(const ftype *__restrict__ beam_dt,
                                   ftype *__restrict__ beam_dE, const int n_rf,
                                   const ftype *__restrict__ voltage, const ftype *__restrict__ omega_RF,
                                   const ftype *__restrict__ phi_RF, const int n_macroparticles,
                                   const ftype acc_kick, const int_vector_t &filter)
{
// KICK
   int k = 0;
   for (int j = 0; j < n_rf; j++) {
      for (const auto &i : filter)
         beam_dE[i] += voltage[k] *
                       vdt::fast_sin(omega_RF[k] * beam_dt[i] + phi_RF[k]);
      k += GP->n_turns;
   }

// SYNCHRONOUS ENERGY CHANGE
   for (const auto &i : filter)
      beam_dE[i] += acc_kick;

}

//drift without periodicity
inline void RingAndRfSection::drift(ftype *__restrict__ beam_dt,
                                    const ftype *__restrict__ beam_dE, const solver_type solver,
                                    const ftype T0, const ftype length_ratio, const int alpha_order,
                                    const ftype eta_zero, const ftype eta_one, const ftype eta_two,
                                    const ftype beta, const ftype energy, const int n_macroparticles)
{

//beam_dt[0] += 0.000001;

   int i;
   ftype T = T0 * length_ratio;

   if (solver == simple) {
      ftype coeff = eta_zero / (beta * beta * energy);
      #pragma omp parallel for
      for (i = 0; i < n_macroparticles; i++)
         beam_dt[i] += T * coeff * beam_dE[i];
   }

   else {
      const ftype coeff = 1. / (beta * beta * energy);
      const ftype eta0 = eta_zero * coeff;
      const ftype eta1 = eta_one * coeff * coeff;
      const ftype eta2 = eta_two * coeff * coeff * coeff;

      if (alpha_order == 1)
         for (i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
      else if (alpha_order == 2)
         for (i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                     - eta1 * beam_dE[i]
                                     * beam_dE[i]) - 1.);
      else
         for (i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                     - eta1 * beam_dE[i]
                                     * beam_dE[i] - eta2
                                     * beam_dE[i]
                                     * beam_dE[i]
                                     * beam_dE[i]) - 1.);
   }

}

// drift with periodicity

inline void RingAndRfSection::drift(ftype *__restrict__ beam_dt,
                                    const ftype *__restrict__ beam_dE, const solver_type solver,
                                    const ftype T0, const ftype length_ratio, const int alpha_order,
                                    const ftype eta_zero, const ftype eta_one, const ftype eta_two,
                                    const ftype beta, const ftype energy, const int n_macroparticles,
                                    const int_vector_t &filter)
{

   ftype T = T0 * length_ratio;

   if (solver == simple) {
      ftype coeff = eta_zero / (beta * beta * energy);

      for (const auto &i : filter)
         beam_dt[i] += T * coeff * beam_dE[i];
   }

   else {
      const ftype coeff = 1. / (beta * beta * energy);
      const ftype eta0 = eta_zero * coeff;
      const ftype eta1 = eta_one * coeff * coeff;
      const ftype eta2 = eta_two * coeff * coeff * coeff;

      if (alpha_order == 1)
         for (const auto &i : filter)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
      else if (alpha_order == 2)
         for (const auto &i : filter)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                     - eta1 * beam_dE[i]
                                     * beam_dE[i]) - 1.);
      else
         for (const auto &i : filter)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                     - eta1 * beam_dE[i]
                                     * beam_dE[i] - eta2
                                     * beam_dE[i]
                                     * beam_dE[i]
                                     * beam_dE[i]) - 1.);
   }

}


void RingAndRfSection::track()
{


   // Determine phase loop correction on RF phase and frequency
   if (PL != NULL && RfP->counter >= PL->delay)
      PL->track();

   if (periodicity) {
      // Change reference of all the particles on the right of the current
      // frame; these particles skip one kick and drift
      //for (int i = 0; i < Beam->n_macroparticles; ++i) {
      if (indices_right_outside.size() > 0) {
         //std::cout << "Found "
         //          << indices_right_outside.size()
         //          << " right outside particles\n";
         for (const auto &i : indices_right_outside)
            Beam->dt[i] -= GP->t_rev[RfP->counter + 1];
      }

      // Synchronize the bunch with the particles that are on the right of
      // the current frame applying kick and drift to the bunch; after that
      // all the particle are in the new updated frame

      kick(indices_inside_frame, RfP->counter);
      drift(indices_inside_frame, RfP->counter + 1);
      //std::cout << "dt[0] : " << Beam->dt[0] << "\n";
      //std::cout << "dE[0] : " << Beam->dE[0] << "\n";
      // find left outside particles and kick, drift them one more time
      //int a = 0;

      indices_left_outside.clear();
      //#pragma omp parallel for reduction(+:a)
      for (int i = 0; i < Beam->n_macroparticles; ++i) {
         if (Beam->dt[i] < 0) {
            indices_left_outside.push_back(i);
         }
      }
      if (indices_left_outside.size() > 0) {
         //std::cout << "Found "
         //          << indices_left_outside.size()
         //          << " left outside particles\n";

         // This will update only the indices_left_outside values
         //  need to test this

         for (const auto &i : indices_left_outside)
            Beam->dt[i] += GP->t_rev[RfP->counter + 1];

         kick(indices_left_outside, RfP->counter);
         drift(indices_left_outside, RfP->counter + 1);
         //std::cout << "dt[0] : " << Beam->dt[0] << "\n";
         //std::cout << "dE[0] : " << Beam->dE[0] << "\n";

      }

      // update inside, right outside particles

      set_periodicity();

   } else {
      kick(RfP->counter);
      drift(RfP->counter + 1);

   }

   if (dE_max > 0)
      horizontal_cut();

   RfP->counter++;
   //std::cout << "insiders : " << indices_inside_frame.size() << "\n";
   //std::cout << "right : " << indices_right_outside.size() << "\n";
   //std::cout << "left : " << indices_left_outside.size() << "\n";
}

inline void RingAndRfSection::horizontal_cut()
{

   for (int i = 0; i < Beam->n_macroparticles; ++i) {
      if (Beam->dE[i] > - dE_max) {
         Beam->dE.erase(Beam->dE.begin() + i);
         Beam->dt.erase(Beam->dt.begin() + i);
         Beam->id.erase(Beam->id.begin() + i);
      }
   }
   Beam->n_macroparticles = Beam->dE.size();

}

RingAndRfSection::RingAndRfSection(solver_type _solver, PhaseLoop *_PhaseLoop,
                                   ftype *_NoiseFB, bool _periodicity, ftype _dE_max,
                                   bool _rf_kick_interp, ftype *_Slices, ftype *_TotalInducedVoltage)
{
   this->elapsed_time = 0;
   this->solver = _solver;
   this->PL = _PhaseLoop;
   this->NoiseFB = _NoiseFB;
   this->periodicity = _periodicity;
   this->dE_max = _dE_max;
   this->rf_kick_interp = _rf_kick_interp;
   this->Slices = _Slices;
   this->TotalInducedVoltage = _TotalInducedVoltage;


   this->acceleration_kick = new ftype[RfP->n_rf * (GP->n_turns)];
   for (int i = 0; i < RfP->n_rf * GP->n_turns; ++i) {
      acceleration_kick[i] = -RfP->E_increment[i];
   }

   if (solver != simple && solver != full) {
      dprintf(
         "ERROR: Choice of longitudinal solver not recognized! Aborting...");
      exit(-1);
   }

   if (GP->alpha_order > 1) {
      solver = full;
   }

   if (periodicity) {
      for (int i = 0; i < Beam->n_macroparticles; i++) {
         if (Beam->dt[i] < 0) {
            dprintf("ERROR: condition Beam.dt >= 0 not true!");
            exit(-1);
         }
      }
      set_periodicity();
      // we will either do this by using vectors or not do it at all if there is no actual need
      GP->t_rev.push_back(GP->t_rev.back());
   }

}

void RingAndRfSection::set_periodicity()
{
   indices_right_outside.clear();
   indices_inside_frame.clear();
   // TODO I dont duplicate the insiders dE, dt
   // as done in the python version
   for (int i = 0; i < Beam->n_macroparticles; ++i) {
      if (Beam->dt[i] > GP->t_rev[RfP->counter + 1]) {
         //std::cout << "Found a right outside particle!\n";
         indices_right_outside.push_back(i);
      } else {
         indices_inside_frame.push_back(i);
      }
   }
}

inline void RingAndRfSection::kick(const int index)
{
   kick(Beam->dt.data(), Beam->dE.data(), RfP->n_rf, &RfP->voltage[index],
        &RfP->omega_RF[index], &RfP->phi_RF[index], Beam->n_macroparticles,
        acceleration_kick[index]);
}

void RingAndRfSection::kick(const int_vector_t &filter, const int index)
{
   kick(Beam->dt.data(), Beam->dE.data(), RfP->n_rf, &RfP->voltage[index],
        &RfP->omega_RF[index], &RfP->phi_RF[index], Beam->n_macroparticles,
        acceleration_kick[index], filter);
}

void RingAndRfSection::drift(const int_vector_t &filter, const int index)
{
   drift(Beam->dt.data(), Beam->dE.data(), solver, GP->t_rev[index], RfP->length_ratio,
         GP->alpha_order, RfP->eta_0(index), RfP->eta_1(index),
         RfP->eta_2(index), RfP->beta(index), RfP->energy(index),
         Beam->n_macroparticles, filter);
}

inline void RingAndRfSection::drift(const int index)
{
   drift(Beam->dt.data(), Beam->dE.data(), solver, GP->t_rev[index], RfP->length_ratio,
         GP->alpha_order, RfP->eta_0(index), RfP->eta_1(index),
         RfP->eta_2(index), RfP->beta(index), RfP->energy(index),
         Beam->n_macroparticles);
}
