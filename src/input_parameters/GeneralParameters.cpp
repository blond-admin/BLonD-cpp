/*
 * GeneralParameters.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include <blond/input_parameters/GeneralParameters.h>

GeneralParameters::GeneralParameters(const uint _n_turns,
                                     f_vector_t &_ring_length,
                                     f_vector_2d_t &_alpha,
                                     const uint _alpha_order,
                                     f_vector_2d_t &_momentum,
                                     const particle_type _particle,
                                     ftype user_mass,
                                     ftype user_charge,
                                     const particle_type _particle2,
                                     ftype user_mass_2,
                                     ftype user_charge_2,
                                     const uint number_of_sections)
{

   this->particle = _particle;
   this->particle_2 = _particle2;
   this->n_sections = number_of_sections;

   if (particle == proton) {
      mass = constant::m_p * constant::c * constant::c / constant::e;
      charge = 1;
   } else if (particle == electron) {
      mass = constant::m_e * constant::c * constant::c / constant::e;
      charge = -1;
   } else if (particle == user_input) {
      mass = user_mass;
      charge = user_charge;
   } else {
      dprintf("ERROR: Particle type not recognized!");
      exit(-1);
   }

   if (particle_2 == none) {
      ;
   } else if (particle_2 == proton) {
      mass2 = constant::m_p * constant::c * constant::c / constant::e;
      charge2 = 1;
   } else if (particle == electron) {
      mass2 = constant::m_e * constant::c * constant::c / constant::e;
      charge2 = -1;
   } else if (particle == user_input) {
      mass2 = user_mass_2;
      charge2 = user_charge_2;
   } else {
      dprintf("ERROR: Second particle type not recognized!");
      exit(-1);
   }

   this->n_turns = _n_turns;
   this->momentum = _momentum;
   this->alpha_order = _alpha_order - 1;
   this->alpha = _alpha;
   this->ring_length = _ring_length;
   this->ring_circumference = std::accumulate(&ring_length[0],
                              &ring_length[n_sections], 0.0);
   this->ring_radius = ring_circumference / (2 * constant::pi);

   if (n_sections > 1) {
      // TODO do some things inside here
      // Should ask danilo about this
      // Danilo told me we could skip this for now
   }

   this->gamma.resize(n_sections, f_vector_t(n_turns + 1));
   this->beta.resize(n_sections, f_vector_t(n_turns + 1));
   this->energy.resize(n_sections, f_vector_t(n_turns + 1));
   this->kin_energy.resize(n_sections, f_vector_t(n_turns + 1));

   const ftype masssq = mass * mass;

   for (uint i = 0; i < n_sections; ++i) {
      for (uint j = 0; j < n_turns + 1; ++j) {
         const ftype momentumsq = momentum[i][j] * momentum[i][j];
         this->beta[i][j] = std::sqrt(1 / (1 + (masssq / momentumsq)));
         this->gamma[i][j] = std::sqrt(1 + (momentumsq / masssq));
         this->energy[i][j] = std::sqrt(masssq + momentumsq);
         this->kin_energy[i][j] = energy[i][j] - mass;
      }
   }

   t_rev.resize(n_turns + 1, 0);

   for (uint i = 0; i < n_sections; ++i)
      for (uint j = 0; j < n_turns + 1; ++j)
         t_rev[j] += ring_length[i] / (beta[i][j] * constant::c);


   cycle_time.resize(n_turns);
   cycle_time[0] = 0;
   for (uint i = 1; i < n_turns; ++i)
      cycle_time[i] = t_rev[i] + cycle_time[i - 1];

   f_rev.resize(n_turns + 1);
   for (uint i = 0; i < n_turns + 1; ++i)
      f_rev[i] = 1 / t_rev[i];

   omega_rev.resize(n_turns + 1);
   for (uint i = 0; i < n_turns + 1; ++i)
      omega_rev[i] = 2 * constant::pi * f_rev[i];

   if (alpha_order > 3) {
      dprintf(
         "WARNING: Momentum compaction factor is implemented only up to 2nd order");
      alpha_order = 3;
   }
   this->eta_0.resize(n_sections, f_vector_t(n_turns + 1));
   this->eta_1.resize(n_sections, f_vector_t(n_turns + 1));
   this->eta_2.resize(n_sections, f_vector_t(n_turns + 1));
   eta_generation();
}


GeneralParameters::~GeneralParameters() {}


void GeneralParameters::eta_generation()
{
   _eta0();
   if (alpha_order > 0)
      _eta1();
   if (alpha_order > 1)
      _eta2();
   if (alpha_order > 2)
      dprintf(
         "WARNING: Momentum compaction factor is implemented only up to 2nd order");
}

void GeneralParameters::_eta0()
{
   //eta_0 = new ftype[n_sections * (n_turns + 1)];
   for (uint i = 0; i < n_sections; ++i)
      for (uint j = 0; j < n_turns + 1; ++j)
         eta_0[i][j] = alpha[i][0] - 1 / (gamma[i][j] * gamma[i][j]);
   //dprintf("eta_0[0] = %lf\n", eta_0[0]);

}

void GeneralParameters::_eta1()
{
   //eta_1 = new ftype[n_sections * (n_turns + 1)];
   for (uint i = 0; i < n_sections; ++i)
      for (uint j = 0; j < n_turns + 1; ++j)
         eta_1[i][j] = 3 * beta[i][j] * beta[i][j]
                       / (2 * gamma[i][j] * gamma[i][j])
                       + alpha[i][1] - alpha[i][0] * eta_0[i][j];


}

void GeneralParameters::_eta2()
{
   //eta_2 = new ftype[n_sections * (n_turns + 1)];
   for (uint i = 0; i < n_sections; ++i)
      for (uint j = 0; j < n_turns + 1; ++j) {
         const ftype betasq = beta[i][j] * beta[i][j];
         ftype gammasq = gamma[i][j] * gamma[i][j];
         eta_1[i][j] = - betasq * (5 * betasq - 1)
                       / (2 * gammasq) + alpha[i][2]
                       - 2 * alpha[i][0] * alpha[i][1]
                       + alpha[i][1] / gammasq
                       + alpha[i][0] * alpha[i][0] * eta_0[i][j]
                       - 3 * betasq * alpha[i][0] / (2 * gammasq);
      }

}



