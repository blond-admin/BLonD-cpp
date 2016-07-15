/*
 * Tracker.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include <blond/trackers/Tracker.h>
#include <blond/math_functions.h>

// Two versions of kick, drift one with periodicity and another without periodiciy
// First go the versions without periodicity

// Kick without periodicity
inline void RingAndRfSection::kick(const ftype *__restrict beam_dt,
                                   ftype *__restrict beam_dE,
                                   const int n_rf,
                                   const ftype *__restrict voltage,
                                   const ftype *__restrict omega_RF,
                                   const ftype *__restrict phi_RF,
                                   const int n_macroparticles,
                                   const ftype acc_kick)
{
   // KICK
   //#pragma omp parallel for collapse(2)
   for (int j = 0; j < n_rf; ++j) {
      #pragma omp parallel for
      for (int i = 0; i < n_macroparticles; ++i) {
         //const ftype a = omega_RF[j] * beam_dt[i] + phi_RF[j];
         beam_dE[i] += voltage[j]
                       * mymath::fast_sin(omega_RF[j] * beam_dt[i]
                                          + phi_RF[j]);
      }
   }

   // SYNCHRONOUS ENERGY CHANGE
   #pragma omp parallel for
   for (int i = 0; i < n_macroparticles; ++i)
      beam_dE[i] += acc_kick;

}


// kick with periodicity

inline void RingAndRfSection::kick(const ftype *__restrict beam_dt,
                                   ftype *__restrict beam_dE,
                                   const int n_rf,
                                   const ftype *__restrict voltage,
                                   const ftype *__restrict omega_RF,
                                   const ftype *__restrict phi_RF,
                                   const int n_macroparticles,
                                   const ftype acc_kick,
                                   const int_vector_t &filter)
{
   // KICK
   //#pragma omp parallel for collapse(2)
   for (int j = 0; j < n_rf; j++) {
      for (const auto &i : filter) {
         const ftype a = omega_RF[j] * beam_dt[i] + phi_RF[j];
         beam_dE[i] += voltage[j] * mymath::fast_sin(a);
      }
   }

   // SYNCHRONOUS ENERGY CHANGE
   for (const auto &i : filter)
      beam_dE[i] += acc_kick;

}

//drift without periodicity
inline void RingAndRfSection::drift(ftype *__restrict beam_dt,
                                    const ftype *__restrict beam_dE,
                                    const solver_type solver,
                                    const ftype T0,
                                    const ftype length_ratio,
                                    const uint alpha_order,
                                    const ftype eta_zero,
                                    const ftype eta_one,
                                    const ftype eta_two,
                                    const ftype beta,
                                    const ftype energy,
                                    const int n_macroparticles)
{

   const ftype T = T0 * length_ratio;

   if (solver == simple) {
      const ftype T_x_coeff = T * eta_zero / (beta * beta * energy);
      #pragma omp parallel for
      for (int i = 0; i < n_macroparticles; i++)
         beam_dt[i] += T_x_coeff * beam_dE[i];
   } else {
      const ftype coeff = 1. / (beta * beta * energy);
      const ftype eta0 = eta_zero * coeff;
      const ftype eta1 = eta_one * coeff * coeff;
      const ftype eta2 = eta_two * coeff * coeff * coeff;

      if (alpha_order == 1)
         #pragma omp parallel for
         for (int i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
      else if (alpha_order == 2)
         #pragma omp parallel for
         for (int i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                     - eta1 * beam_dE[i]
                                     * beam_dE[i]) - 1.);
      else
         #pragma omp parallel for
         for (int i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                     - eta1 * beam_dE[i]
                                     * beam_dE[i] - eta2
                                     * beam_dE[i]
                                     * beam_dE[i]
                                     * beam_dE[i]) - 1.);
   }

}

// drift with periodicity

inline void RingAndRfSection::drift(ftype *__restrict beam_dt,
                                    const ftype *__restrict beam_dE,
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
	auto GP = Context::GP;
	auto RfP = Context::RfP;
	auto Beam = Context::Beam;


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
      for (uint i = 0; i < Beam->n_macroparticles; ++i) {
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
	auto Beam = Context::Beam;

   for (uint i = 0; i < Beam->n_macroparticles; ++i) {
      if (Beam->dE[i] > - dE_max) {
         Beam->dE.erase(Beam->dE.begin() + i);
         Beam->dt.erase(Beam->dt.begin() + i);
         Beam->id.erase(Beam->id.begin() + i);
      }
   }
   Beam->n_macroparticles = Beam->dE.size();

}

RingAndRfSection::RingAndRfSection(solver_type _solver,
                                   PhaseLoop *_PhaseLoop,
                                   LHCNoiseFB *_NoiseFB,
                                   bool _periodicity,
                                   ftype _dE_max,
                                   bool _rf_kick_interp,
                                   Slices *_Slices,
                                   TotalInducedVoltage *_TotalInducedVoltage)
{
	auto GP = Context::GP;
	auto RfP = Context::RfP;
	auto Beam = Context::Beam;

	this->elapsed_time = 0;
   this->solver = _solver;
   this->PL = _PhaseLoop;
   this->noiseFB = _NoiseFB;
   this->periodicity = _periodicity;
   this->dE_max = _dE_max;
   this->rf_kick_interp = _rf_kick_interp;
   this->slices = _Slices;
   this->totalInducedVoltage = _TotalInducedVoltage;


   this->acceleration_kick.resize(GP->n_turns + 1);
   // = new ftype[RfP->n_rf * (GP->n_turns)];
   for (uint i = 0; i < RfP->E_increment.size(); ++i) {
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
      for (uint i = 0; i < Beam->n_macroparticles; i++) {
         if (Beam->dt[i] < 0) {
            dprintf("ERROR: condition Beam.dt >= 0 not true!");
            exit(-1);
         }
      }
      set_periodicity();
      GP->t_rev.push_back(GP->t_rev.back());
   }

   // vol = new ftype[RfP->n_rf];
   // omeg = new ftype[RfP->n_rf];
   // phi = new ftype[RfP->n_rf];

}

RingAndRfSection::~RingAndRfSection()
{
   // delete[] vol;
   // delete[] omeg;
   // delete[] phi;
}


void RingAndRfSection::set_periodicity()
{
	auto GP = Context::GP;
	auto RfP = Context::RfP;
	auto Beam = Context::Beam;

   indices_right_outside.clear();
   indices_inside_frame.clear();
   // TODO I am not duplicating the insiders dE, dt
   // as done in the python version
   for (uint i = 0; i < Beam->n_macroparticles; ++i) {
      if (Beam->dt[i] > GP->t_rev[RfP->counter + 1]) {
         indices_right_outside.push_back(i);
      } else {
         indices_inside_frame.push_back(i);
      }
   }
}

inline void RingAndRfSection::kick(const uint index)
{
	auto RfP = Context::RfP;
	auto Beam = Context::Beam;

   auto vol = new ftype[RfP->n_rf];
   auto omeg = new ftype[RfP->n_rf];
   auto phi = new ftype[RfP->n_rf];

   for (uint i = 0; i < RfP->n_rf; ++i) {
      vol[i] = RfP->voltage[i][index];
      omeg[i] = RfP->omega_RF[i][index];
      phi[i] = RfP->phi_RF[i][index];
   }

   kick(Beam->dt.data(),
        Beam->dE.data(),
        RfP->n_rf,
        vol,
        omeg,
        phi,
        Beam->n_macroparticles,
        acceleration_kick[index]);

   delete[] vol;
   delete[] omeg;
   delete[] phi;
}

void RingAndRfSection::kick(const int_vector_t &filter, const uint index)
{
	auto RfP = Context::RfP;
	auto Beam = Context::Beam;

   auto vol = new ftype[RfP->n_rf];
   auto omeg = new ftype[RfP->n_rf];
   auto phi = new ftype[RfP->n_rf];

   for (uint i = 0; i < RfP->n_rf; ++i) {
      vol[i] = RfP->voltage[i][index];
      omeg[i] = RfP->omega_RF[i][index];
      phi[i] = RfP->phi_RF[i][index];
   }

   kick(Beam->dt.data(),
        Beam->dE.data(),
        RfP->n_rf,
        vol,
        omeg,
        phi,
        Beam->n_macroparticles,
        acceleration_kick[index],
        filter);

   delete[] vol;
   delete[] omeg;
   delete[] phi;
}

void RingAndRfSection::drift(const int_vector_t &filter, const uint index)
{
	auto GP = Context::GP;
	auto RfP = Context::RfP;
	auto Beam = Context::Beam;

   drift(Beam->dt.data(),
         Beam->dE.data(),
         solver,
         GP->t_rev[index],
         RfP->length_ratio,
         GP->alpha_order,
         RfP->eta_0(index),
         RfP->eta_1(index),
         RfP->eta_2(index),
         RfP->beta(index),
         RfP->energy(index),
         Beam->n_macroparticles,
         filter);
}

inline void RingAndRfSection::drift(const uint index)
{
	auto GP = Context::GP;
	auto RfP = Context::RfP;
	auto Beam = Context::Beam;

   drift(Beam->dt.data(),
         Beam->dE.data(), solver,
         GP->t_rev[index],
         RfP->length_ratio,
         GP->alpha_order,
         RfP->eta_0(index),
         RfP->eta_1(index),
         RfP->eta_2(index),
         RfP->beta(index),
         RfP->energy(index),
         Beam->n_macroparticles);
}
