/*
 * PhaseLoop.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: kiliakis
 */

#include "PhaseLoop.h"

PhaseLoop::PhaseLoop(f_vector_t PL_gain,
                     ftype window_coefficient,
                     uint _delay,
                     PhaseNoise *phaseNoise,
                     LHCNoiseFB *LHCNoiseFB)
{
}


void PhaseLoop::beam_phase()
{
   /*
    *Beam phase measured at the main RF frequency and phase. The beam is
    convolved with the window function of the band-pass filter of the
    machine. The coefficients of sine and cosine components determine the
    beam phase, projected to the range -Pi/2 to 3/2 Pi. Note that this beam
    phase is already w.r.t. the instantaneous RF phase.*
    */

   // Main RF frequency at the present turn
   //omega_RF = self.rf_params.omega_RF[0,self.rf_params.counter[0]]
   ftype omega_RF = RfP->omega_RF[RfP->idx][RfP->counter];
   ftype phi_RF = RfP->phi_RF[RfP->idx][RfP->counter];
   // Convolve with window function
   //
   ftype *base = new ftype[Slice->n_slices];
   ftype *array = new ftype[Slice->n_slices];

   #pragma omp parallel for
   for (int i = 0; i < (int)Slice->n_slices; ++i) {
      const ftype a = alpha * Slice->bin_centers[i];
      base[i] = std::exp(a) * Slice->n_macroparticles[i];
   }

   #pragma omp parallel for
   for (int i = 0; i < (int)Slice->n_slices; ++i) {
      const ftype a = omega_RF * Slice->bin_centers[i] + phi_RF;
      array[i] = base[i] * std::sin(a);
   }
   ftype scoeff = mymath::trapezoid(array, Slice->bin_centers.data(),
                                    Slice->n_slices);

   #pragma omp parallel for
   for (int i = 0; i < (int)Slice->n_slices; ++i) {
      const ftype a = omega_RF * Slice->bin_centers[i] + phi_RF;
      array[i] = base[i] * std::cos(a);
   }

   ftype ccoeff = mymath::trapezoid(array, Slice->bin_centers.data(),
                                    Slice->n_slices);

   phi_beam = std::atan(scoeff / ccoeff) + constant::pi;

   delete[] base;
   delete[] array;

}


void PhaseLoop::phase_difference()
{
   /*
    Phase difference between beam and RF phase of the main RF system.
    Optional: add RF phase noise through dphi directly.*
    */

   // Correct for design stable phase
   uint counter = RfP->counter;
   dphi = phi_beam - RfP->phi_s[counter];
   // Possibility to add RF phase noise through the PL
   if (RFnoise != NULL) {
      if (noiseFB != NULL) {
         dphi += noiseFB->fX * RFnoise->fDphi[counter];
      } else {
         dphi += RFnoise->fDphi[counter];
      }
   }
}


void PhaseLoop::radial_steering_from_freq()
{
   // Frequency and phase change for the current turn due to the
   // radial steering program.

   const uint counter = RfP->counter;

   const auto radial_steering_domega_RF =
      - RfP->omega_RF_d[0][counter] * RfP->eta_0(counter)
      / GP->alpha[0][0] * reference / GP->ring_radius;

   for (uint i = 0; i < RfP->n_rf; ++i) {
      RfP->omega_RF[i][counter] +=
         radial_steering_domega_RF * RfP->harmonic[i][counter]
         / RfP->harmonic[0][counter];
   }

   // Update the RF phase of all systems for the next turn
   // Accumulated phase offset due to PL in each RF system
   for (uint i = 0; i < RfP->n_rf; ++i) {
      RfP->dphi_RF_steering[i] += 2 * constant::pi * RfP->harmonic[i][counter]
                                  * (RfP->omega_RF[i][counter]
                                     - RfP->omega_RF_d[i][counter])
                                  / RfP->omega_RF_d[i][counter];
   }

   // Total phase offset
   for (uint i = 0; i < RfP->n_rf; ++i) {
      RfP->phi_RF[i][counter] += RfP->dphi_RF_steering[i];
   }

}


void PhaseLoop::radial_difference()
{
   // Radial difference between beam and design orbit.*
   uint counter = RfP->counter;
   uint n = 0;
   ftype sum = 0;
   //ftype array[GP->n_turns];
   for (uint i = 0; i < GP->n_turns; ++i) {
      if (Beam->dt[i] > Slice->bin_centers.front()
            and Beam->dt[i] < Slice->bin_centers.back()) {
         sum += Beam->dE[i];
         n++;
      }
   }
   auto average_dE = n > 0 ? sum / n : 0.0;
   // std::cout << "average_dE : " << average_dE << "\n";
   drho = GP->alpha[0][0] * GP->ring_radius * average_dE
          / (GP->beta[0][counter] * GP->beta[0][counter]
             * GP->energy[0][counter]);
   // std::cout << "drho : " << drho << "\n";

}


void PhaseLoop::default_track()
{
   /*
    Calculate PL correction on main RF frequency depending on machine.
    Update the RF phase and frequency of the next turn for all systems.
    */

   uint counter = RfP->counter + 1;
   //uint turns = GP->n_turns;
   // Update the RF frequency of all systems for the next turn
   for (uint i = 0; i < RfP->n_rf; ++i) {
      RfP->omega_RF[i][counter] += domega_RF * RfP->harmonic[i][counter]
                                   / RfP->harmonic[0][counter];
   }

   //Update the RF phase of all systems for the next turn
   //Accumulated phase offset due to PL in each RF system
   for (uint i = 0; i < RfP->n_rf; ++i) {
      RfP->dphi_RF[i] += 2 * constant::pi * RfP->harmonic[i][counter]
                         * (RfP->omega_RF[i][counter]
                            - RfP->omega_RF_d[i][counter])
                         / RfP->omega_RF_d[i][counter];
   }

   // Total phase offset
   for (uint i = 0; i < RfP->n_rf; ++i) {
      RfP->phi_RF[i][counter] += RfP->dphi_RF[i];
   }

}



LHC::LHC(f_vector_t PL_gain,
         ftype SL_gain,
         ftype window_coefficient,
         PhaseNoise *phaseNoise ,
         LHCNoiseFB *LHCNoiseFB,
         uint _delay)
{

   // General Initializations
   this->delay = _delay;
   this->alpha = window_coefficient;
   this->gain = PL_gain;
   this->RFnoise = phaseNoise;
   this->noiseFB = LHCNoiseFB;
   // End of general initializations

   this->gain2 = SL_gain;
   this->lhc_y = 0;
   this->lhc_a.resize(GP->n_turns + 1, 0);
   this->lhc_t.resize(GP->n_turns + 1, 0);

   if (gain2 != 0) {
      for (uint i = 0; i < GP->n_turns + 1; ++i) {
         lhc_a[i] = 5.25 - RfP->omega_s0[i] / (constant::pi * 40);
      }
      for (uint i = 0; i < GP->n_turns + 1; ++i) {
         lhc_t[i] = (2 * constant::pi * RfP->Qs[i] * sqrt(lhc_a[i]))
                    / sqrt(1 + gain[i] / gain2
                           * sqrt((1 + 1 / lhc_a[i])
                                  / (1 + lhc_a[i])));
      }
   }
}


LHC::~LHC() {}


void LHC::track()
{
   /*
    Calculation of the LHC RF frequency correction from the phase difference
    between beam and RF (actual synchronous phase). The transfer function is

    .. math::
    \\Delta \\omega_{rf}^{PL} = - g_{PL} (\\Delta\\varphi_{PL} + \\phi_{N})

    where the phase noise for the controlled blow-up can be optionally
    activated.
    Using 'gain2', a synchro loop can be activated in addition to remove
    long-term frequency drifts:

    .. math::
    \\Delta \\omega_{rf}^{SL} = - g_{SL} (y + a \\Delta\\varphi_{rf}) ,

    where we use the recursion

    .. math::
    y_{n+1} = (1 - \\tau) y_n + (1 - a) \\tau \\Delta\\varphi_{rf} ,

    with a and \tau being defined through the synchrotron frequency f_s and
    the synchrotron tune Q_s as

    .. math::
    a (f_s) \\equiv 5.25 - \\frac{f_s}{\\pi 40~\\text{Hz}} ,

    .. math::
    \\tau(f_s) \\equiv 2 \\pi Q_s \\sqrt{ \\frac{a}{1 + \\frac{g_{PL}}{g_{SL}} \\sqrt{\\frac{1 + 1/a}{1 + a}} }}
    */
   uint counter = RfP->counter;
   ftype dphi_RF = RfP->dphi_RF[0];

   beam_phase();
   phase_difference();

   // Frequency correction from phase loop and synchro loop

   domega_RF = -gain[counter] * dphi - gain2
               * (lhc_y + lhc_a[counter] * (dphi_RF + reference));

   // Update recursion variable
   lhc_y = (1 - lhc_t[counter]) * lhc_y + (1 - lhc_a[counter])
           * lhc_t[counter] * (dphi_RF + reference);

   default_track();
}


PSB::PSB(f_vector_t PL_gain,
         f_vector_t _RL_gain,
         ftype _PL_period,
         ftype _RL_period,
         f_vector_t coefficients,
         ftype window_coefficient,
         PhaseNoise *phaseNoise,
         LHCNoiseFB *LHCNoiseFB,
         uint _delay)
{

   // General Initializations
   this->delay = _delay;
   this->alpha = window_coefficient;
   this->gain = PL_gain;
   this->RFnoise = phaseNoise;
   this->noiseFB = LHCNoiseFB;
   // End of general initializations


   // figure out what to do with these pairs
   this->gain2.resize(2);// = new ftype[2];
   this->dt.resize(2);// = new uint[2];

   if (_RL_gain.empty()) {
      gain2[0] = 0;
      gain2[1] = 0;
   } else {
      gain2[0] = _RL_gain[0];
      gain2[1] = _RL_gain[1];
   }

   if (_PL_period == 0)
      dt[0] = 10e-6;
   else
      dt[0] = _PL_period;

   if (_RL_period == 0)
      dt[1] = 7;
   else
      dt[1] = _RL_period;

   this->PL_counter = 1;
   this->on_time.push_back(0);

   precalculate_time();

   // TODO coefficients
   // How many can i have?

   //*Memory of previous phase correction, for phase loop.*
   dphi_av = 0;
   dphi_av_prev = 0;
   //*Memory of previous relative radial correction, for rad loop.*
   drho_prev = 0;
   //*Accumulated time for radial loop*
   t_accum = 0;
   //*Phase loop frequency correction [1/s]*
   domega_PL = 0;
   //*Radial loop frequency correction [1/s]*
   domega_RL = 0;
   //domega_RF = 0;
}


PSB::~PSB() {}


void PSB::track()
{
   /*
    Calculation of the PSB RF frequency correction from the phase difference
    between beam and RF (actual synchronous phase). The transfer function is

    .. math::
    \\Delta \\omega_{RF} = g(t) \\frac{a_0 \\Delta\\Phi_{PL}^2 + a_1 \\Delta\\Phi_{PL} + a_2 }{b_0 \\Delta\\Phi_{PL}^2 + b_1 \\Delta\\Phi_{PL} + b_2}

    Input g through gain and [a_0, a_1, a_2, b_0, b_1, b_2] through coefficients.
    */

   // Average phase error while frequency is updated
   uint counter = RfP->counter;

   beam_phase();
   phase_difference();
   dphi_av += dphi;
   t_accum += GP->t_rev[counter];

   //dprintf("Before if counter is %d\n", counter);
   // Phase loop active on certain turns

   if (counter == on_time[PL_counter] && counter > delay) {
      dphi_av /= (on_time[PL_counter]
                  - on_time[PL_counter - 1]);

      domega_PL = 0.998 * domega_PL
                  - gain[counter]
                  * (dphi_av - dphi_av_prev + reference);

      // Update averaging variables
      dphi_av_prev = dphi_av;
      dphi_av = 0;

      // Add correction from radial loop
      if (PL_counter % static_cast<int>(dt[1]) == 0) {

         drho = (RfP->omega_RF[0][counter] - RfP->omega_RF_d[0][counter])
                / (RfP->omega_RF_d[0][counter]
                   * (1 / (GP->alpha[0][0] * RfP->gamma(counter)
                           * RfP->gamma(counter)) - 1))
                + reference;

         domega_RL = domega_RL - gain2[0] * (drho - drho_prev)
                     - gain2[1] * drho * t_accum;

         drho_prev = drho;
         t_accum = 0;
      }

      // Counter to pick the next time step when the PL will be active
      PL_counter++;
   }
   // Apply frequency correction
   domega_RF = domega_PL + domega_RL;

   default_track();

}


void PSB::precalculate_time()
{
   /*
    For machines like the PSB, where the PL acts only in certain time
    intervals, pre-calculate on which turns to act.
    */
   uint n = delay + 1;

   while (n < GP->t_rev.size()) {
      auto summa = 0.0;
      while (summa < dt[0]) {
         if (n < GP->t_rev.size()) {
            summa += GP->t_rev[n];
            n++;
         } else {
            on_time.push_back(0);
            return;
         }
      }
      on_time.push_back(n - 1);
   }

}


LHC_F::LHC_F(ftype PL_gain,
             ftype window_coefficient,
             ftype FL_gain,
             PhaseNoise *phaseNoise,
             LHCNoiseFB *LHCNoiseFB,
             uint delay)
{
   this->gain = PL_gain;
   this->alpha = window_coefficient;
   this->delay = delay;
   this->RFnoise = phaseNoise;
   this->noiseFB = LHCNoiseFB;
   this->gain2 = FL_gain;
   // domega_RF.resize(gain.size());
}


LHC_F::~LHC_F() {}


void LHC_F::track()
{
   uint counter = RfP->counter;

   beam_phase();
   phase_difference();

   // Frequency correction from phase loop and frequency loop
   const auto factor = RfP->omega_RF[0][counter]
                       - RfP->omega_RF_d[0][counter]
                       + reference;
   // for (int i = 0; i < (int)gain.size(); ++i) {
   domega_RF = - gain * dphi
               - gain2 * factor;
   // }

   default_track();
}


SPS_RL::SPS_RL(ftype PL_gain,
               ftype window_coefficient,
               ftype RL_gain,
               PhaseNoise *phaseNoise,
               LHCNoiseFB *LHCNoiseFB,
               uint delay)
{
   this->gain = PL_gain;
   this->alpha = window_coefficient;
   this->delay = delay;
   this->RFnoise = phaseNoise;
   this->noiseFB = LHCNoiseFB;
   this->gain2 = RL_gain;
}


SPS_RL::~SPS_RL() {}


// TODO Test this function
void SPS_RL::track()
{
   uint counter = RfP->counter;

   if (reference != 0)
      radial_steering_from_freq();

   beam_phase();
   phase_difference();
   radial_difference();

   // Frequency correction from phase loop and radial loop
   const auto factor = mymath::sign(RfP->eta_0(counter))
                       * (reference - drho) / GP->ring_radius;

   // for (int i = 0; i < (int)gain.size(); ++i) {
   domega_RF = - gain * dphi
               - gain2 * factor;
   // }


   default_track();
}
