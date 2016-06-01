/*
 * GeneralParameters.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */
#include <blond/input_parameters/GeneralParameters.h>

#include <numeric>
#include <blond/utilities.h>
#include <blond/constants.h>
using namespace blond;

GeneralParameters::~GeneralParameters()
{
   util::delete_array(alpha);
   //util::delete_array (momentum);
   util::delete_array(ring_length);
   util::delete_array(beta);
   util::delete_array(gamma);
   util::delete_array(energy);
   util::delete_array(kin_energy);
   util::delete_array(cycle_time);
   util::delete_array(f_rev);
   util::delete_array(omega_rev);
   util::delete_array(eta_0);
   util::delete_array(eta_1);
   util::delete_array(eta_2);
   //delete t_rev;
   
}


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
   for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
      int j = i / (n_turns + 1);
      eta_0[i] = alpha[j] - 1 / (gamma[i] * gamma[i]);
   }
   //dprintf("eta_0[0] = %lf\n", eta_0[0]);

}

void GeneralParameters::_eta1()
{
   //eta_1 = new ftype[n_sections * (n_turns + 1)];
   for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
      int j = i / (n_turns + 1);
      eta_1[i] = 3 * beta[i] * beta[i] / (2 * gamma[i] * gamma[i])
                 + alpha[j + 1] - alpha[j] * eta_0[i];
   }

}

void GeneralParameters::_eta2()
{
   //eta_2 = new ftype[n_sections * (n_turns + 1)];
   for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
      int j = i / (n_turns + 1);
      ftype betasq = beta[i] * beta[i];
      ftype gammasq = gamma[i] * gamma[i];
      eta_1[i] = -betasq * (5 * betasq - 1) / (2 * gammasq) + alpha[j + 2]
                 - 2 * alpha[j] * alpha[j + 1] + alpha[j + 1] / gammasq
                 + alpha[j] * alpha[j] * eta_0[i]
                 - 3 * betasq * alpha[j] / (2 * gammasq);
   }
}

GeneralParameters::GeneralParameters(const int _n_turns, ftype *_ring_length,
                                     ftype *_alpha, const int _alpha_order, ftype *_momentum,
                                     const particle_type _particle, ftype user_mass, ftype user_charge,
                                     particle_type _particle2, ftype user_mass_2, ftype user_charge_2,
                                     int number_of_sections)
{

   //global++;
   //dprintf("Global variable is %d\n", global);
   this->n_turns = _n_turns;
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

   // We don't have tuple momentum == we don't need cumulative_times
   // assuming that momentum is a 2d array

   this->momentum = _momentum;

   this->alpha_order = _alpha_order;

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

   this->gamma = new ftype[n_sections * (n_turns + 1)];
   this->beta = new ftype[n_sections * (n_turns + 1)];
   this->energy = new ftype[n_sections * (n_turns + 1)];
   this->kin_energy = new ftype[n_sections * (n_turns + 1)];

   ftype masssq = mass * mass;

   for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
      ftype momentumsq = momentum[i] * momentum[i];
      this->beta[i] = sqrt(1 / (1 + (masssq / momentumsq)));
      this->gamma[i] = sqrt(1 + (momentumsq / masssq));
      this->energy[i] = sqrt(masssq + momentumsq);
      this->kin_energy[i] = energy[i] - mass;
   }

   for (int i = 0; i < n_turns + 1; ++i)
      t_rev.push_back(0);
   //dump(ring_length, 1, "ring_length\n");

   for (int j = 0; j < n_turns + 1; ++j)
      for (int k = 0; k < n_sections; ++k)
         t_rev[j] += ring_length[k] / (beta[k * (n_turns + 1) + j] * constant::c);

   //dump(t_rev, 10, "t_rev\n");

   this->cycle_time = new ftype[n_turns];
   cycle_time[0] = 0;//t_rev[0];
   for (int i = 1; i < n_turns; ++i)
      cycle_time[i] = t_rev[i] + cycle_time[i - 1];

   //dump(cycle_time, 10, "cycle_time\n");

   this->f_rev = new ftype[n_turns + 1];
   for (int i = 0; i < n_turns + 1; ++i)
      f_rev[i] = 1 / t_rev[i];

   this->omega_rev = new ftype[n_turns + 1];
   for (int i = 0; i < n_turns + 1; ++i)
      omega_rev[i] = 2 * constant::pi * f_rev[i];

   if (alpha_order > 3) {
      dprintf(
         "WARNING: Momentum compaction factor is implemented only up to 2nd order");
      alpha_order = 3;
   }
   this->eta_0 = new ftype[n_sections * (n_turns + 1)];
   this->eta_1 = new ftype[n_sections * (n_turns + 1)];
   this->eta_2 = new ftype[n_sections * (n_turns + 1)];

   eta_generation();
}

