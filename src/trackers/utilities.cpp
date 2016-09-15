/*
 * utilities.cpp
 *
 *  Created on: Sep 5, 2016
 *      Author: kiliakis
 */

#include <cmath>
#include <iostream>
#include <blond/trackers/utilities.h>
#include <blond/constants.h>
#include <blond/math_functions.h>


std::vector<bool> is_in_separatrix(const GeneralParameters *GP,
                                   const RfParameters *RfP,
                                   const Beams *Beam,
                                   const f_vector_t &dt,
                                   const f_vector_t &dE,
                                   f_vector_t total_voltage)
{
    /*
    Condition for being inside the separatrix.
    For the time being, for single RF section only or from total voltage.
    Single RF sinusoidal.
    Uses beta, energy averaged over the turn.
    To be generalized.
    */

    if (GP->n_sections > 1)
        std::cerr << "WARNING: is_in_separatrix is not yet properly computed"
                  << "for several sections!\n";
    if (RfP->n_rf > 1)
        std::cerr << "WARNING: is_in_separatrix will be calculated for"
                  << "the first harmonic only!\n";

    const int counter = RfP->counter;
    const ftype dt_sep = (constant::pi - RfP->phi_s[counter]
                          - RfP->phi_RF[0][counter])
                         / RfP->omega_RF[0][counter];

    const ftype Hsep = std::abs(hamiltonian(GP, RfP, Beam, dt_sep, 0.0));

    auto temp = hamiltonian(GP, RfP, Beam, dt.data(),
                            dE.data(), dt.size(),
                            total_voltage);

    std::vector<bool> isin(temp.size());
    const int size = isin.size();
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        isin[i] = std::abs(temp[i]) < Hsep;

    return isin;
}



// TODO fix the issue with eta tracking with dE vector and alpha_order > 1
f_vector_t hamiltonian(const GeneralParameters *GP,
                       const RfParameters *RfP,
                       const Beams *Beam,
                       const ftype *__restrict dt,
                       const ftype *__restrict dE,
                       const int size,
                       const f_vector_t total_voltage)
{
    if (GP->n_sections > 1)
        std::cerr << "WARNING: The Hamiltonian is not yet properly computed"
                  << "for several sections!\n";
    if (RfP->n_rf > 1)
        std::cerr << "WARNING: The Hamiltonian will be calculated for"
                  << "the first harmonic only!\n";

    const int counter = RfP->counter;
    const ftype h0 = RfP->harmonic[0][counter];
    ftype v0;

    if (total_voltage.empty())
        v0 = RfP->voltage[0][counter];
    else
        v0 = total_voltage[counter];
    v0 *= GP->charge;

    // TODO it should be dE instead of 0.0
    const ftype c1 = RfP->eta_tracking(Beam, counter, 0.0)
                     * constant::pi * constant::c
                     / (GP->ring_circumference * RfP->beta(counter) *
                        RfP->energy(counter));

    const ftype c2 = constant::c * RfP->beta(counter) * v0
                     / (h0 * GP->ring_circumference);
    const ftype phi_s = RfP->phi_s[counter];

    f_vector_t phi_b(size);
    const ftype omega_rf0 = RfP->omega_RF[0][counter];
    const ftype phi_rf0 = RfP->phi_RF[0][counter];

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        phi_b[i] = omega_rf0 * dt[i] + phi_rf0;

    const ftype eta0 = RfP->eta_0(counter);

    if (eta0 < 0)
        #pragma omp parallel for
        for (int i = 0; i < size; i++)
            phi_b[i] = phase_modulo_below_transition(phi_b[i]);
    else if (eta0 > 0)
        #pragma omp parallel for
        for (int i = 0; i < size; i++)
            phi_b[i] = phase_modulo_above_transition(phi_b[i]);

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        phi_b[i] =  c1 * dE[i] * dE[i]
                    + c2 * (mymath::fast_cos(phi_b[i]) - mymath::fast_cos(phi_s)
                            + (phi_b[i] - phi_s) * mymath::fast_sin(phi_s));


    return phi_b;
}


ftype hamiltonian(const GeneralParameters *GP,
                  const RfParameters *RfP,
                  const Beams *Beam,
                  const ftype dt,
                  const ftype dE,
                  const f_vector_t total_voltage)
{
    if (GP->n_sections > 1)
        std::cerr << "WARNING: The Hamiltonian is not yet properly computed"
                  << "for several sections!\n";
    if (RfP->n_rf > 1)
        std::cerr << "WARNING: The Hamiltonian will be calculated for"
                  << "the first harmonic only!\n";

    const int counter = RfP->counter;
    const ftype h0 = RfP->harmonic[0][counter];
    ftype v0;

    if (total_voltage.empty())
        v0 = RfP->voltage[0][counter];
    else
        v0 = total_voltage[counter];
    v0 *= GP->charge;

    const ftype c1 = RfP->eta_tracking(Beam, counter, dE)
                     * constant::pi * constant::c
                     / (GP->ring_circumference * RfP->beta(counter) *
                        RfP->energy(counter));

    const ftype c2 = constant::c * RfP->beta(counter) * v0
                     / (h0 * GP->ring_circumference);

    const ftype phi_s = RfP->phi_s[counter];
    ftype phi_b = RfP->omega_RF[0][counter] * dt + RfP->phi_RF[0][counter];
    const ftype eta0 = RfP->eta_0(counter);

    if (eta0 < 0)
        phi_b = phase_modulo_below_transition(phi_b);
    else if (eta0 > 0)
        phi_b = phase_modulo_above_transition(phi_b);


    return c1 * dE * dE
           + c2 * (mymath::fast_cos(phi_b) - mymath::fast_cos(phi_s)
                   + (phi_b - phi_s) * mymath::fast_sin(phi_s));
}