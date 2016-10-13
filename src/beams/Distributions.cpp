/*
 * Beams.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */
#include <blond/beams/Distributions.h>
#include <blond/configuration.h>
#include <blond/constants.h>
#include <blond/utilities.h>
#include <blond/globals.h>
#include <iostream>
#include <blond/trackers/Tracker.h>
#include <map>

using namespace std;

void matched_from_line_density(FullRingAndRf *full_ring,
                               map<string, string> line_density_opt,
                               string main_harmonic,
                               string plot,
                               string figdir,
                               string half_option,
                               map<string, string> extraVoltageDict,
                               int n_iterations_input,
                               int seed
                              )
{

}


void matched_from_distribution_density(FullRingAndRf *full_ring,
                                       map<string, string> distribution_opt,
                                       string main_harmonic,
                                       int n_iterations_input,
                                       map<string, string> extraVoltageDict,
                                       int seed
                                      )
{

}


f_vector_t distribution_density_function(const f_vector_t &action_array,
        const string &dist_type, const double length, double exponent)
{
    f_vector_t ret(action_array.size());
    if (dist_type == "binomial") {
        for (uint i = 0; i < action_array.size(); i++) {
            if (action_array[i] > length)
                ret[i] = 0.;
            else
                ret[i] = pow((1. - action_array[i] / length), exponent);
        }
        return ret;
    } else if (dist_type == "waterbag") {
        for (uint i = 0; i < action_array.size(); i++) {
            if (action_array[i] > length)
                ret[i] = 0.;
            else
                ret[i] = 1.;
        }
        return ret;
    } else if (dist_type == "parabolic_amplitude") {
        for (uint i = 0; i < action_array.size(); i++) {
            if (action_array[i] > length)
                ret[i] = 0.;
            else
                ret[i] = (1 - action_array[i] / length);
        }
        return ret;
    } else if (dist_type == "parabolic_line") {
        exponent = 0.5;
        for (uint i = 0; i < action_array.size(); i++) {
            if (action_array[i] > length)
                ret[i] = 0.;
            else
                ret[i] = pow((1. - action_array[i] / length), exponent);
        }
        return ret;
    } else if (dist_type == "gaussian") {
        for (uint i = 0; i < action_array.size(); i++)
            ret[i] = exp(-2. * action_array[i] / length);
    } else {
        cerr << "[distribution_density_function] The dist_type was not recognized\n";
        exit(-1);
    }

    return ret;

}



f_vector_t line_density_function(const f_vector_t &coord_array,
                                 const string &dist_type,
                                 const double bunch_length,
                                 const double bunch_position,
                                 double exponent)
{
    f_vector_t ret(coord_array.size());
    if (dist_type == "binomial") {
        for (uint i = 0; i < coord_array.size(); i++) {
            if (abs(coord_array[i] - bunch_position) > bunch_length / 2.)
                ret[i] = 0.;
            else
                ret[i] = pow(1 -
                             pow((coord_array[i] - bunch_position)
                                 / (bunch_length / 2.), 2.),
                             exponent + 0.5);
        }
        return ret;
    } else if (dist_type == "waterbag") {
        exponent = 0.;
        for (uint i = 0; i < coord_array.size(); i++) {
            if (abs(coord_array[i] - bunch_position) > bunch_length / 2.)
                ret[i] = 0.;
            else
                ret[i] = pow(1 -
                             pow((coord_array[i] - bunch_position)
                                 / (bunch_length / 2.), 2.),
                             exponent + 0.5);
        }
        return ret;
    } else if (dist_type == "parabolic_amplitude") {
        exponent = 1.;
        for (uint i = 0; i < coord_array.size(); i++) {
            if (abs(coord_array[i] - bunch_position) > bunch_length / 2.)
                ret[i] = 0.;
            else
                ret[i] = pow(1 -
                             pow((coord_array[i] - bunch_position)
                                 / (bunch_length / 2.), 2.),
                             exponent + 0.5);
        }
        return ret;
    } else if (dist_type == "parabolic_line") {
        exponent = 0.5;
        for (uint i = 0; i < coord_array.size(); i++) {
            if (abs(coord_array[i] - bunch_position) > bunch_length / 2.)
                ret[i] = 0.;
            else
                ret[i] = pow(1 -
                             pow((coord_array[i] - bunch_position)
                                 / (bunch_length / 2.), 2.),
                             exponent + 0.5);
        }
        return ret;
    } else if (dist_type == "gaussian") {
        const double sigma = bunch_length / 4.;
        for (uint i = 0; i < coord_array.size(); i++)
            ret[i] = exp(- pow(coord_array[i] - bunch_position, 2)
                         / (2 * sigma * sigma));
    } else if (dist_type == "cosine_squared") {
        for (uint i = 0; i < coord_array.size(); i++) {
            if (abs(coord_array[i] - bunch_position) > bunch_length / 2.)
                ret[i] = 0.;
            else
                ret[i] = pow(cos(constant::pi
                                 * (coord_array[i] - bunch_position)
                                 / bunch_length),
                             2);
        }
        return ret;
    } else {
        cerr << "[line_density_function] The dist_type was not recognized\n";
        exit(-1);
    }

    return ret;
}



void longitudinal_bigaussian(ftype sigma_dt, ftype sigma_dE,
                             int seed, bool reinsertion)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    if (GP->n_sections > 1) {
        cerr << "WARNING: longitudinal_bigaussian is not yet properly\n"
             << "computed for several sections!\n";
    }
    if (RfP->n_rf > 1) {
        cerr << "WARNING: longitudinal_bigaussian for multiple RF\n"
             << "is not yet implemented\n";
    }

    uint counter = RfP->counter;
    ftype harmonic = RfP->harmonic[0][counter];
    ftype energy = GP->energy[0][counter];
    ftype beta = GP->beta[0][counter];
    ftype omega_RF = RfP->omega_RF[0][counter];
    ftype phi_s = RfP->phi_s[counter];
    ftype phi_RF = RfP->phi_RF[0][counter];
    ftype eta0 = RfP->eta_0(counter);

    if (sigma_dE == 0) {
        auto voltage = GP->charge * RfP->voltage[0][counter];
        auto phi_b = omega_RF * sigma_dt + phi_s;
        sigma_dE =
            sqrt(voltage * energy * beta * beta *
                 (cos(phi_b) - cos(phi_s) + (phi_b - phi_s) * sin(phi_s)) /
                 (constant::pi * harmonic * eta0));
    }

    Beam->sigma_dE = sigma_dE;
    Beam->sigma_dt = sigma_dt;

    if (seed < 0) {
        f_vector_t random;
        util::read_vector_from_file(random, TEST_FILES "/normal_distribution.dat");
        for (uint i = 0; i < Beam->n_macroparticles; ++i) {
            ftype r = random[i % random.size()];
            if (eta0 > 0)
                Beam->dt[i] = sigma_dt * r
                              + (phi_s - phi_RF) / omega_RF;
            else
                Beam->dt[i] = sigma_dt * r
                              + (phi_s - phi_RF - constant::pi) / omega_RF;
            Beam->dE[i] = sigma_dE * r;
        }
    } else {
        default_random_engine generator(seed);
        normal_distribution<ftype> distribution(0.0, 1.0);
        for (uint i = 0; i < Beam->n_macroparticles; ++i) {
            ftype r = distribution(generator);
            if (eta0 > 0)
                Beam->dt[i] = sigma_dt * r + (phi_s - phi_RF) / omega_RF;
            else
                Beam->dt[i] =
                    sigma_dt * r + (phi_s - phi_RF - constant::pi) / omega_RF;
            r = distribution(generator);
            Beam->dE[i] = sigma_dE * r;
        }
    }

    // TODO if reinsertion == true
    if (reinsertion) {
        ;
    }
}