/*
 * GeneralParameters.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include <blond/constants.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <numeric>
#include <blond/vector_math.h>
#include <blond/utilities.h>

GeneralParameters::GeneralParameters(
    const int _n_turns, f_vector_t &_ring_length, f_vector_2d_t &_alpha,
    const int _alpha_order, f_vector_2d_t &_momentum,
    const particle_t _particle, ftype user_mass, ftype user_charge,
    const particle_t _particle2, ftype user_mass_2, ftype user_charge_2,
    const int number_of_sections)
{

    particle = _particle;
    particle_2 = _particle2;
    n_sections = number_of_sections;

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
        std::cerr << "ERROR: Particle type not recognized!\n";
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
        std::cerr << "ERROR: Second particle type not recognized!\n";
        exit(-1);
    }

    n_turns = _n_turns;
    momentum = _momentum;
    alpha_order = _alpha_order - 1;
    alpha = _alpha;
    ring_length = _ring_length;
    ring_circumference = std::accumulate(ALL(ring_length), 0.0);
    ring_radius = ring_circumference / (2 * constant::pi);

    if (n_sections > 1) {
        // TODO do some things inside here
        // Should ask danilo about this
        // Danilo told me we could skip this for now
    }

    gamma.resize(n_sections, f_vector_t(n_turns + 1));
    beta.resize(n_sections, f_vector_t(n_turns + 1));
    energy.resize(n_sections, f_vector_t(n_turns + 1));
    kin_energy.resize(n_sections, f_vector_t(n_turns + 1));

    const ftype masssq = mass * mass;

    for (int i = 0; i < n_sections; ++i) {
        for (int j = 0; j < n_turns + 1; ++j) {
            const ftype momentumsq = momentum[i][j] * momentum[i][j];
            beta[i][j] = std::sqrt(1 / (1 + (masssq / momentumsq)));
            gamma[i][j] = std::sqrt(1 + (momentumsq / masssq));
            energy[i][j] = std::sqrt(masssq + momentumsq);
            kin_energy[i][j] = energy[i][j] - mass;
        }
    }

    t_rev.resize(n_turns + 1, 0);

    for (int i = 0; i < n_sections; ++i)
        t_rev += ring_length[i] / constant::c / beta[i];

    cycle_time.resize(n_turns);
    cycle_time[0] = 0;
    for (int i = 1; i < n_turns; ++i)
        cycle_time[i] = t_rev[i] + cycle_time[i - 1];

    f_rev = 1.0 / t_rev;

    omega_rev = 2. * constant::pi * f_rev;

    if (alpha_order > 3) {
        dprintf(
            "WARNING: Momentum compaction factor is implemented only up to 2nd "
            "order");
        alpha_order = 3;
    }
    eta_0.resize(n_sections, f_vector_t(n_turns + 1));
    eta_1.resize(n_sections, f_vector_t(n_turns + 1));
    eta_2.resize(n_sections, f_vector_t(n_turns + 1));
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
            "WARNING: Momentum compaction factor is implemented only up to 2nd "
            "order");
}

void GeneralParameters::_eta0()
{
    for (int i = 0; i < n_sections; ++i)
        for (int j = 0; j < n_turns + 1; ++j)
            eta_0[i][j] = alpha[i][0] - 1 / (gamma[i][j] * gamma[i][j]);
}

void GeneralParameters::_eta1()
{
    for (int i = 0; i < n_sections; ++i)
        for (int j = 0; j < n_turns + 1; ++j)
            eta_1[i][j] =
                3 * beta[i][j] * beta[i][j] / (2 * gamma[i][j] * gamma[i][j]) +
                alpha[i][1] - alpha[i][0] * eta_0[i][j];
}

void GeneralParameters::_eta2()
{
    for (int i = 0; i < n_sections; ++i)
        for (int j = 0; j < n_turns + 1; ++j) {
            const ftype betasq = beta[i][j] * beta[i][j];
            ftype gammasq = gamma[i][j] * gamma[i][j];
            eta_1[i][j] = -betasq * (5 * betasq - 1) / (2 * gammasq) +
                          alpha[i][2] - 2 * alpha[i][0] * alpha[i][1] +
                          alpha[i][1] / gammasq +
                          alpha[i][0] * alpha[i][0] * eta_0[i][j] -
                          3 * betasq * alpha[i][0] / (2 * gammasq);
        }
}
