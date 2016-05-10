/*
 * RfParameters.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include "RfParameters.h"

RfParameters::~RfParameters()
{
   util::delete_array(this->E_increment);
   util::delete_array(this->phi_s);
   util::delete_array(this->Qs);
   util::delete_array(this->omega_RF_d);
   util::delete_array(this->omega_s0);
   util::delete_array(this->phi_RF);
   util::delete_array(this->dphi_RF);
   util::delete_array(this->dphi_RF_steering);
   util::delete_array(this->t_RF);
   util::delete_array(this->omega_RF);
   util::delete_array(this->harmonic);
   util::delete_array(this->voltage);
   util::delete_array(this->phi_offset);
   util::delete_array(this->phi_noise);
}

/*
 :How to use RF programs:

 - For 1 RF system and constant values of V, h or phi, just input the single value
 - For 1 RF system and varying values of V, h or phi, input an array of n_turns values
 - For several RF systems and constant values of V, h or phi, input lists of single values
 - For several RF systems and varying values of V, h or phi, input lists of arrays of n_turns values
 */
// RfParameters == RfSectionParameters
// completely removed accelerating_systems
RfParameters::RfParameters(int _n_rf, ftype *_harmonic, ftype *_voltage,
                           ftype *_phi_offset, ftype *_phi_noise, ftype *_omega_RF,
                           int _section_index, accelerating_systems_type accelerating_systems)
{

   this->counter = 0;
   this->section_index = _section_index - 1;

   this->n_rf = _n_rf;
   this->harmonic = _harmonic;
   this->voltage = _voltage;
   this->phi_offset = _phi_offset;
   this->phi_noise = NULL;

   this->section_length = GP->ring_length[section_index];
   this->length_ratio = section_length / GP->ring_circumference;

   // wiped out all the imports
   this->E_increment = new ftype[GP->n_turns];
   // Don't have n_turns +1 cause of the way np.diff works
   for (int j = 0; j < GP->n_turns; ++j) {
      E_increment[j] = GP->energy[j + 1] - GP->energy[j];
   }

   this->phi_s = new ftype[GP->n_turns + 1];
   calc_phi_s(accelerating_systems);

   this->Qs = new ftype[(GP->n_turns + 1)];
   for (int i = 0; i < GP->n_turns + 1; ++i) {
      Qs[i] = sqrt(
                 harmonic[i] * GP->charge * voltage[i]
                 * fabs(eta_0(i) * cos(phi_s[i]))
                 / (2 * constant::pi * GP->beta[i] * GP->beta[i] * GP->energy[i]));
      //if (i < 3)
      //dprintf("%lf \n", harmonic[i] * GP->charge * voltage[i]);
      //dprintf("%.4e \n", std::abs(eta_0(i) * cos(phi_s[i])));

   }
   //dump(Qs, 10, "Qs\n");

   this->omega_s0 = new ftype[(GP->n_turns + 1)];
   for (int i = 0; i < (GP->n_turns + 1); ++i) {
      this->omega_s0[i] = Qs[i] * GP->omega_rev[i];
   }
   //dump(omega_s0, 10, "omega_s0\n");


   this->omega_RF_d = new ftype[n_rf * (GP->n_turns + 1)];
   for (int i = 0; i < n_rf * (GP->n_turns + 1); ++i) {
      this->omega_RF_d[i] = 2 * constant::pi * GP->beta[i] * constant::c * harmonic[i]
                            / GP->ring_circumference;
   }
   //dump(omega_RF_d, 10, "omega_RF_d\n");
   //dprintf("ring_circumference %.12lf\n", GP->ring_circumference);
   //dprintf("pi %.12lf\n", pi);
   //dprintf("c %.12lf\n", c);
   //dprintf("pi*c %.12lf\n", pi*c);
   /*
    this->omega_RF = new ftype[n_rf * (GP->n_turns + 1)];
    if (omega_RF == NULL) {
    dprintf("It was null!\n");
    for (int i = 0; i < n_rf * (GP->n_turns + 1); ++i) {
    omega_RF[i] = omega_RF_d[i];
    }
    }
    */
   if (_omega_RF == NULL) {
      this->omega_RF = new ftype[n_rf * (GP->n_turns + 1)];
      std::copy(&omega_RF_d[0], &omega_RF_d[n_rf * (GP->n_turns + 1)],
                omega_RF);
   } else {
      this->omega_RF = _omega_RF;
   }

   this->phi_RF = new ftype[n_rf * (GP->n_turns + 1)];
   //this->phi_RF = (ftype *) aligned_malloc(
   //    sizeof(ftype) * n_rf * (GP->n_turns + 1));
   std::copy(&phi_offset[0], &phi_offset[n_rf * (GP->n_turns + 1)], phi_RF);

   this->dphi_RF = new ftype[n_rf];
   std::fill_n(dphi_RF, n_rf, 0);

   this->dphi_RF_steering = new ftype[n_rf];
   std::fill_n(dphi_RF_steering, n_rf, 0);

   this->t_RF = new ftype[GP->n_turns + 1];

   for (int i = 0; i < GP->n_turns + 1; ++i) {
      t_RF[i] = 2 * constant::pi / omega_RF[i];
   }
}

ftype RfParameters::eta_tracking(const Beams *beam, const int counter,
                                 const ftype dE)
{
   ftype eta = 0;
   if (GP->alpha_order == 1)
      eta = eta_0(counter);
   else {
      ftype delta = dE / ((GP->beta[0]) * (GP->beta[0]) * GP->energy[0]);
      eta += eta_0(counter) * 1;
      if (GP->alpha_order > 0)
         eta += eta_1(counter) * delta;
      if (GP->alpha_order > 1)
         eta += eta_2(counter) * delta * delta;
      if (GP->alpha_order > 2)
         dprintf(
            "WARNING: Momentum compaction factor is implemented only up to 2nd order");
   }
   return eta;

}

ftype RfParameters::eta_0(const int i)
{
   return GP->eta_0[section_index * (GP->n_turns + 1) + i];
}

ftype RfParameters::eta_1(const int i)
{
   return GP->eta_1[section_index * (GP->n_turns + 1) + i];
}

ftype RfParameters::eta_2(const int i)
{
   return GP->eta_2[section_index * (GP->n_turns + 1) + i];
}

ftype RfParameters::beta(const int i)
{
   return GP->beta[section_index * (GP->n_turns + 1) + i];
}

ftype RfParameters::gamma(const int i)
{
   return GP->gamma[section_index * (GP->n_turns + 1) + i];

}

ftype RfParameters::energy(const int i)
{
   return GP->energy[section_index * (GP->n_turns + 1) + i];

}

ftype RfParameters::momentum(const int i)
{
   return GP->momentum[section_index * (GP->n_turns + 1) + i];

}

int RfParameters::sign_eta_0(const int i)
{
   if (GP->eta_0[section_index * (GP->n_turns + 1) + i] > 0)
      return 1;
   else if (GP->eta_0[section_index * (GP->n_turns + 1) + i] == 0)
      return 0;
   else
      return -1;
}

// TODO do I need to pass RfParameters object?
// Do I call this function from somewhere else?
void RfParameters::calc_phi_s(accelerating_systems_type acc_sys)
{
   /*
    | *The synchronous phase calculated from the rate of momentum change.*
    | *Below transition, for decelerating bucket: phi_s is in (-Pi/2,0)*
    | *Below transition, for accelerating bucket: phi_s is in (0,Pi/2)*
    | *Above transition, for accelerating bucket: phi_s is in (Pi/2,Pi)*
    | *Above transition, for decelerating bucket: phi_s is in (Pi,3Pi/2)*
    | *The synchronous phase is calculated at a certain moment.*
    | *Uses beta, energy averaged over the turn.*
    */
   int n_turns = GP->n_turns;
   //ftype eta0 = rf_params->eta0;
   if (acc_sys == as_single) {

      ftype *denergy = new ftype[n_turns + 1];
      for (int j = 0; j < n_turns; ++j)
         denergy[j] = E_increment[j];
      denergy[n_turns] = E_increment[n_turns - 1];

      ftype *acceleration_ratio = new ftype[n_turns + 1];
      for (int i = 0; i < n_turns + 1; ++i)
         acceleration_ratio[i] = denergy[i] / (GP->charge * voltage[0 + i]);

      for (int i = 0; i < n_turns + 1; ++i)
         if (acceleration_ratio[i] > 1 || acceleration_ratio[i] < -1)
            dprintf(
               "Warning!!! Acceleration is not possible (momentum increment is too big or voltage too low) at index %d\n",
               i);

      for (int i = 0; i < n_turns + 1; ++i)
         phi_s[i] = asin(acceleration_ratio[i]);

      //ftype *eta0_middle_points = new ftype[n_turns +1];
      for (int i = 0; i < n_turns; ++i) {
         ftype middle = (eta_0(i) + eta_0(i + 1)) / 2;
         if (middle > 0)
            phi_s[i] = constant::pi - phi_s[i];
         else
            phi_s[i] = constant::pi + phi_s[i];
      }
      if (eta_0(n_turns) > 0)
         phi_s[n_turns] = constant::pi - phi_s[n_turns];
      else
         phi_s[n_turns] = constant::pi + phi_s[n_turns];

      delete[] denergy;
      delete[] acceleration_ratio;

      return;

   } else if (acc_sys == all) {
      /*
       In this case, all the RF systems are accelerating, phi_s is
       calculated accordingly with respect to the fundamental frequency (the minimum
       of the potential well is taken)
       */

      ftype *transition_phase_offset = new ftype[n_turns + 1];
      for (int i = 0; i < n_turns + 1; ++i) {
         phi_s[i] = 0;
         if (eta_0(i) > 0)
            transition_phase_offset[i] = constant::pi;
         else
            transition_phase_offset[i] = 0;
      }
      ftype phase_array[1000];
      mymath::linspace(phase_array, -constant::pi * 1.2, constant::pi * 1.2, 1000);

      for (int i = 0; i < n_turns; ++i) {
         ftype totalRF[1000] = { };
         for (int j = 0; j < n_rf; j++) {
            int row = j * (n_turns + 1);
            ftype min = harmonic[mymath::min(harmonic, n_rf, n_turns + 1)];

            for (int k = 0; k < 1000; ++k) {
               totalRF[k] += voltage[row + i + 1]
                             * vdt::fast_sin(
                                (harmonic[row + i + 1] / min)
                                * (phase_array[k]
                                   + transition_phase_offset[i
                                         + 1])
                                + phi_offset[row + i + 1]);
            }
         }
         int transition_factor = transition_phase_offset[i] == 0 ? +1 : -1;
         //dump(totalRF, 10, "totalRF\n");

         ftype potential_well[1000] = { 0 };
         ftype *f = new ftype[1000];
         for (int k = 0; k < 1000; ++k) {
            f[k] = totalRF[k] - E_increment[i] / GP->charge;
         }

         //dump(f, 10, "f\n");

         //dprintf("dx %.12lf\n", phase_array[1] - phase_array[0]);

         ftype *trap = mymath::cum_trapezoid(f,
                                             phase_array[1] - phase_array[0], 1000);

         for (int k = 0; k < 1000; ++k) {
            potential_well[k] = transition_factor * trap[k];
         }
         //dump(potential_well, 10, "potential_well\n");

         // TODO why mean here? line BLonD-minimal::rf_parameter.py:334
         phi_s[i + 1] = phase_array[mymath::min(potential_well, 1000, 1)]
                        + transition_phase_offset[i + 1];

      }
      delete[] transition_phase_offset; 
      phi_s[0] = phi_s[1];
      //dump(phi_s, 10, "phi_s\n");

      return;

   } else if (acc_sys == first) {
      /*
       Only the first RF system is accelerating, so we have to correct the
       phi_offset of the other rf_systems such that p_increment relates
       only to the first RF
       */
      ;
   } else {
      dprintf(
         "Did not recognize the option accelerating_systems in calc_phi_s function\n");
      exit(-1);
   }

   // TODO, how can we get here?
   // TODO why n_turns and not n_turns+1?

   if (eta_0(0) > 0) {
      for (int i = 0; i < n_turns; ++i)
         phi_s[i] = constant::pi;
   } else {
      for (int i = 0; i < n_turns; ++i)
         phi_s[i] = 0;
   }

}

