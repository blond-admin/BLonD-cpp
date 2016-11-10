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

using namespace std;


vector<int> is_in_separatrix(const GeneralParameters *GP,
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
        cerr << "WARNING: is_in_separatrix is not yet properly computed"
             << "for several sections!\n";
    if (RfP->n_rf > 1)
        cerr << "WARNING: is_in_separatrix will be calculated for"
             << "the first harmonic only!\n";

    const int counter = RfP->counter;
    const double dt_sep = (constant::pi - RfP->phi_s[counter]
                           - RfP->phi_rf[0][counter])
                          / RfP->omega_rf[0][counter];

    const double Hsep = abs(hamiltonian(GP, RfP, Beam, dt_sep, 0.0));

    auto temp = hamiltonian(GP, RfP, Beam, dt.data(),
                            dE.data(), dt.size(),
                            total_voltage);

    // vector<bool> is not thread safe!!
    vector<int> isin(temp.size());
    const int size = isin.size();

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        isin[i] = abs(temp[i]) < Hsep;

    return isin;
}



// TODO fix the issue with eta tracking with dE vector and alpha_order > 1
f_vector_t hamiltonian(const GeneralParameters *GP,
                       const RfParameters *RfP,
                       const Beams *Beam,
                       const double *__restrict dt,
                       const double *__restrict dE,
                       const int size,
                       const f_vector_t total_voltage)
{
    if (GP->n_sections > 1)
        cerr << "WARNING: The Hamiltonian is not yet properly computed"
             << "for several sections!\n";
    if (RfP->n_rf > 1)
        cerr << "WARNING: The Hamiltonian will be calculated for"
             << "the first harmonic only!\n";

    const int counter = RfP->counter;
    const double h0 = RfP->harmonic[0][counter];
    double v0;

    if (total_voltage.empty())
        v0 = RfP->voltage[0][counter] * GP->charge;
    else
        v0 = total_voltage[counter] * GP->charge;

    // TODO it should be dE instead of 0.0
    const double c1 = RfP->eta_tracking(Beam, counter, 0.0)
                      * constant::pi * constant::c
                      / (GP->ring_circumference * RfP->beta[counter] *
                         RfP->energy[counter]);

    const double c2 = constant::c * RfP->beta[counter] * v0
                      / (h0 * GP->ring_circumference);
    const double phi_s = RfP->phi_s[counter];

    const double omega_rf0 = RfP->omega_rf[0][counter];
    const double phi_rf0 = RfP->phi_rf[0][counter];

    f_vector_t phi_b(size, phi_rf0);

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        phi_b[i] += omega_rf0 * dt[i];

    const double eta0 = RfP->eta_0[counter];

    if (eta0 < 0)
        #pragma omp parallel for
        for (int i = 0; i < size; i++)
            phi_b[i] = phase_modulo_below_transition(phi_b[i]);
    else if (eta0 > 0)
        #pragma omp parallel for
        for (int i = 0; i < size; i++)
            phi_b[i] = phase_modulo_above_transition(phi_b[i]);

    const double cos_phi_s = mymath::fast_cos(phi_s);
    const double sin_phi_s = mymath::fast_sin(phi_s);

    // auto f1 = [ = ](double phib, double de) {
    //     return c1 * de * de + c2
    //            * (mymath::fast_cos(phib) - cos_phi_s + (phib - phi_s) * sin_phi_s);
    // };
    // transform(phi_b.begin(), phi_b.end(), &dE[0], phi_b.begin(), f1);

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        phi_b[i] = c1 * dE[i] * dE[i]
                   + c2 * (mymath::fast_cos(phi_b[i]) - cos_phi_s
                           + (phi_b[i] - phi_s) * sin_phi_s);


    return phi_b;
}


double hamiltonian(const GeneralParameters *GP,
                   const RfParameters *RfP,
                   const Beams *Beam,
                   const double dt,
                   const double dE,
                   const f_vector_t total_voltage)
{
    if (GP->n_sections > 1)
        cerr << "WARNING: The Hamiltonian is not yet properly computed"
             << "for several sections!\n";
    if (RfP->n_rf > 1)
        cerr << "WARNING: The Hamiltonian will be calculated for"
             << "the first harmonic only!\n";

    const int counter = RfP->counter;
    const double h0 = RfP->harmonic[0][counter];
    double v0;

    if (total_voltage.empty())
        v0 = RfP->voltage[0][counter];
    else
        v0 = total_voltage[counter];
    v0 *= GP->charge;

    const double c1 = RfP->eta_tracking(Beam, counter, dE)
                      * constant::pi * constant::c
                      / (GP->ring_circumference * RfP->beta[counter] *
                         RfP->energy[counter]);

    const double c2 = constant::c * RfP->beta[counter] * v0
                      / (h0 * GP->ring_circumference);

    const double phi_s = RfP->phi_s[counter];
    double phi_b = RfP->omega_rf[0][counter] * dt + RfP->phi_rf[0][counter];
    const double eta0 = RfP->eta_0[counter];

    if (eta0 < 0)
        phi_b = phase_modulo_below_transition(phi_b);
    else if (eta0 > 0)
        phi_b = phase_modulo_above_transition(phi_b);


    return c1 * dE * dE
           + c2 * (mymath::fast_cos(phi_b) - mymath::fast_cos(phi_s)
                   + (phi_b - phi_s) * mymath::fast_sin(phi_s));
}


void minmax_location(const f_vector_t &x, const f_vector_t &f,
                     f_vector_t &min_x_position, f_vector_t &max_x_position,
                     f_vector_t &min_values, f_vector_t &max_values)
{
    min_x_position.clear(); max_x_position.clear();
    min_values.clear(); max_values.clear();
    
    f_vector_t f_derivative;
    adjacent_difference(f.begin(), f.end(), back_inserter(f_derivative));
    f_derivative.erase(f_derivative.begin());

    f_vector_t x_derivative;
    for (uint i = 0; i < x.size() - 1; i++)
        x_derivative.push_back(x[i] + (x[1] - x[0]) / 2);


    f_derivative = mymath::interp(x, x_derivative, f_derivative);

    f_vector_t f_derivative_second;
    adjacent_difference(f_derivative.begin(), f_derivative.end(),
                        back_inserter(f_derivative_second));
    f_derivative_second.erase(f_derivative_second.begin());


    f_derivative_second = mymath::interp(x, x_derivative, f_derivative_second);

    f_vector_t f_derivative_zeros;
    for (uint i = 0; i < f_derivative.size() - 1; i++) {
        if (f_derivative[i] == 0. || (f_derivative[i + 1] / f_derivative[i] < 0.))
            f_derivative_zeros.push_back(i);
    }
    for (const auto &i : f_derivative_zeros)
        if (f_derivative_second[i] > 0)
            min_x_position.push_back((x[i + 1] + x[i]) / 2.);
        else
            max_x_position.push_back((x[i + 1] + x[i]) / 2.);

    min_values = mymath::interp(min_x_position, x, f);

    max_values = mymath::interp(max_x_position, x, f);

}

void potential_well_cut(const f_vector_t &theta_coord_array,
                        const f_vector_t &potential_array,
                        f_vector_t &theta_coord_sep,
                        f_vector_t &potential_well_sep)
{
    theta_coord_sep.clear();
    potential_well_sep.clear();

    f_vector_t min_theta_positions, max_theta_positions;
    f_vector_t min_potential_values, max_potential_values;

    minmax_location(theta_coord_array, potential_array,
                    min_theta_positions, max_theta_positions,
                    min_potential_values, max_potential_values);

    int n_minima = min_theta_positions.size();
    int n_maxima = max_theta_positions.size();
    // cout << "minima: " << n_minima << "\n";
    // cout << "maxima: " << n_maxima << "\n";
    if (n_minima == 0) { // tested
        cerr << "[potential_well_cut] The potential well has no minima...\n";
        exit(-1);
    }
    if (n_minima > n_maxima && n_maxima == 1) {
        cerr << "[potential_well_cut] The potential well has more minima,\n"
             << " and only one maximum\n";
        exit(-1);
    }
    if (n_maxima == 0) {
        cout << "[potential_well_cut] Warning: The maximum of the potential \n"
             << "well could not be found... You may reconsider the options to \n"
             << "calculate the potential well as the main harmonic is probably \n"
             << "not the expected one. You may also increase the percentage of \n"
             << "margin to compute the potential well.\n"
             << "The full potential well will be taken\n";
    } else if (n_maxima == 1) { // tested
        if (min_theta_positions[0] > max_theta_positions[0]) {
            for (uint i = 0; i < potential_array.size(); ++i) {
                if (potential_array[i] < max_potential_values[0]
                        && theta_coord_array[i] > max_theta_positions[0]) {
                    theta_coord_sep.push_back(theta_coord_array[i]);
                    potential_well_sep.push_back(potential_array[i]);
                }
            }
            if (potential_array.back() < potential_array[0]) {
                cerr << "[potential_well_cut] The potential well is not well defined.\n"
                     << "You may reconsider the options to calculate \n"
                     << "the potential well as the main harmonic is \n"
                     << "probably not the expected one.\n";
                exit(-1);
            }
        } else {
            for (uint i = 0; i < potential_array.size(); ++i) {
                if (potential_array[i] < max_potential_values[0]
                        && theta_coord_array[i] < max_theta_positions[0]) {
                    theta_coord_sep.push_back(theta_coord_array[i]);
                    potential_well_sep.push_back(potential_array[i]);
                }
            }
            if (potential_array.back() < potential_array[0]) {
                cerr << "[potential_well_cut] The potential well is not well defined.\n"
                     << "You may reconsider the options to calculate \n"
                     << "the potential well as the main harmonic is \n"
                     << "probably not the expected one.\n";
                exit(-1);
            }
        }
    } else if (n_maxima == 2) { // tested
        auto lower_maximum_value = *min_element(ALL(max_potential_values));
        auto higher_maximum_value = *max_element(ALL(max_potential_values));
        f_vector_t lower_maximum_theta, higher_maximum_theta;
        for (uint i = 0; i < max_theta_positions.size(); i++) {
            if (max_potential_values[i] == lower_maximum_value)
                lower_maximum_theta.push_back(max_theta_positions[i]);
            if (max_potential_values[i] == higher_maximum_value)
                higher_maximum_theta.push_back(max_theta_positions[i]);
        }

        if (lower_maximum_theta.size() == 2) {
            for (uint i = 0; i < potential_array.size(); ++i) {
                if (potential_array[i] < lower_maximum_value
                        && theta_coord_array[i] > lower_maximum_theta[0]
                        && theta_coord_array[i] < lower_maximum_theta[1]) {
                    theta_coord_sep.push_back(theta_coord_array[i]);
                    potential_well_sep.push_back(potential_array[i]);
                }
            }
        } else if (min_theta_positions[0] > lower_maximum_theta[0]) {
            for (uint i = 0; i < potential_array.size(); ++i) {
                if (potential_array[i] < lower_maximum_value
                        && theta_coord_array[i] > lower_maximum_theta[0]
                        && theta_coord_array[i] < higher_maximum_theta[0]) {
                    theta_coord_sep.push_back(theta_coord_array[i]);
                    potential_well_sep.push_back(potential_array[i]);
                }
            }
        } else {
            for (uint i = 0; i < potential_array.size(); ++i) {
                if (potential_array[i] < lower_maximum_value
                        && theta_coord_array[i] < lower_maximum_theta[0]
                        && theta_coord_array[i] > higher_maximum_theta[0]) {
                    theta_coord_sep.push_back(theta_coord_array[i]);
                    potential_well_sep.push_back(potential_array[i]);
                }
            }
        }

    } else { // tested
        auto left_max_theta = *min_element(max_theta_positions.begin(),
                                           max_theta_positions.end());
        auto right_max_theta = *max_element(max_theta_positions.begin(),
                                            max_theta_positions.end());

        f_vector_t left_max_value, right_max_value;
        for (uint i = 0; i < max_theta_positions.size(); i++) {
            if (max_theta_positions[i] == left_max_theta)
                left_max_value.push_back(max_potential_values[i]);
            if (max_theta_positions[i] == right_max_theta)
                right_max_value.push_back(max_potential_values[i]);
        }
        auto separatrix_value = min(
                                    *min_element(left_max_value.begin(), left_max_value.end()),
                                    *min_element(right_max_value.begin(), right_max_value.end()));

        for (uint i = 0; i < theta_coord_array.size(); ++i) {
            if (theta_coord_array[i] > left_max_theta
                    && theta_coord_array[i] < right_max_theta
                    && potential_array[i] < separatrix_value) {
                theta_coord_sep.push_back(theta_coord_array[i]);
                potential_well_sep.push_back(potential_array[i]);
            }
        }

    }

}