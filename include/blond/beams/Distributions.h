/*
 * Distributions.h
 *
 *  Created on: Mar 16, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_DISTRIBUTIONS_H_
#define BEAMS_DISTRIBUTIONS_H_

#include <blond/configuration.h>
#include <blond/constants.h>
#include <blond/globals.h>
#include <blond/python.h>
#include <blond/trackers/Tracker.h>
#include <map>
#include <cmath>
#include <random>
#include <stdlib.h>

inline void longitudinal_bigaussian(ftype sigma_dt, ftype sigma_dE = 0,
                                    int seed = 0, bool reinsertion = false)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    if (GP->n_sections > 1) {
        dprintf(
            "WARNING: longitudinal_bigaussian is not yet properly computed for "
            "several sections!");
    }
    if (RfP->n_rf > 1) {
        dprintf("WARNING: longitudinal_bigaussian for multiple RF is not yet "
                "implemented");
    }

    uint counter = RfP->counter;
    ftype harmonic = RfP->harmonic[0][counter];
    ftype energy = GP->energy[0][counter];
    ftype beta = GP->beta[0][counter];
    ftype omega_RF = RfP->omega_RF[0][counter];
    ftype phi_s = RfP->phi_s[counter];
    ftype phi_RF = RfP->phi_RF[0][counter];

    ftype voltage, eta0 = 0.0, phi_b;
    if (sigma_dE == 0) {
        voltage = GP->charge * RfP->voltage[0][counter];
        eta0 = RfP->eta_0(counter);
        phi_b = omega_RF * sigma_dt + phi_s;
        sigma_dE =
            sqrt(voltage * energy * beta * beta *
                 (cos(phi_b) - cos(phi_s) + (phi_b - phi_s) * sin(phi_s)) /
                 (constant::pi * harmonic * eta0));
    }

    Beam->sigma_dE = sigma_dE;
    Beam->sigma_dt = sigma_dt;
    // std::cout << sigma_dE << "\n";
    // std::cout << sigma_dt << "\n";
    // std::cout << (phi_s - phi_RF) / omega_RF << "\n";
    if (seed < 0) {
        for (uint i = 0; i < Beam->n_macroparticles; ++i) {
            ftype r = 1.0 * (i + 1) / Beam->n_macroparticles;
            // ftype r = distribution(generator);
            Beam->dt[i] = sigma_dt * r + (phi_s - phi_RF) / omega_RF;
            // r = 1.0 * rand() / RAND_MAX;
            // r = distribution(generator);
            Beam->dE[i] = sigma_dE * r;
            // dprintf("beam_dE: %.8lf \n", Beam->dE[i]);
        }
    } else {
        std::default_random_engine generator(seed);
        std::normal_distribution<ftype> distribution(0.0, 1.0);
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
        // for (uint i = 0; i < Beam->n_macroparticles; ++i)
        // Beam->dE[i] = sigma_dE * distribution(generator);
    }

    // TODO if reinsertion == true
    if (reinsertion) {
        ;
    }
}




void matched_from_line_density(ftype beta, ftype energy, ftype charge,
                               int n_macroparticles, f_vector_t &dt, f_vector_t &dE,
                               ftype eta_0, ftype t_rev_0,
                               f_vector_t &potential_well_array,
                               f_vector_t &time_coord_array,
                               std::map<std::string, std::string> line_density_opt,
                               std::string main_harmonic = "lowest_freq",
                               std::string figdir = "fig",
                               std::string half_option = "first",
                               std::map<std::string, std::string> extraVoltageDict =
                                   std::map<std::string, std::string>(),
                               int n_iterations_input = 100,
                               int seed = 0
                              )
{
    python::initialize();

    auto pFunc = python::import("distributions", "matched_from_line_density");

    auto pBeta = python::convert_double(beta);
    auto pEnergy = python::convert_double(energy);
    auto pCharge = python::convert_double(charge);
    auto pNMacroparticles = python::convert_int(n_macroparticles);
    auto pEta0 = python::convert_double(eta_0);
    auto pTRev0 = python::convert_double(t_rev_0);
    auto pMainHarmonic = python::convert_string(main_harmonic);
    auto pFigDir = python::convert_string(figdir);
    auto pHalfOption = python::convert_string(half_option);
    auto pNIterationsInput = python::convert_int(n_iterations_input);
    auto pSeed = python::convert_int(seed);
    auto pDT = python::convert_double_array(dt.data(), dt.size());
    auto pDE = python::convert_double_array(dE.data(), dE.size());

    auto pPotentialWell = python::convert_double_array(potential_well_array.data(),
                          potential_well_array.size());

    auto pTimeCoord = python::convert_double_array(time_coord_array.data(),
                      time_coord_array.size());

    auto pLineDensityOpt = python::convert_dictionary(line_density_opt);
    auto pExtraVoltageDict = python::convert_dictionary(extraVoltageDict);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pBeta, pEnergy,
                                            pCharge, pNMacroparticles,
                                            pDT, pDE, pEta0, pTRev0,
                                            pPotentialWell, pTimeCoord,
                                            pLineDensityOpt, pMainHarmonic,
                                            pFigDir, pHalfOption, pExtraVoltageDict,
                                            pNIterationsInput, pSeed, NULL);
    assert(ret);


    python::finalize();
}


void matched_from_line_density(FullRingAndRf *full_ring,
                               std::map<std::string, std::string> line_density_opt =
                                   std::map<std::string, std::string>(),
                               std::string main_harmonic = "lowest_freq",
                               std::string figdir = "fig",
                               std::string half_option = "first",
                               std::map<std::string, std::string> extraVoltageDict =
                                   std::map<std::string, std::string>(),
                               int n_iterations_input = 100,
                               int seed = 0
                              )
{
    auto GP = Context::GP;
    auto Beam = Context::Beam;

    int n_points_potential = int(1e4);

    full_ring->potential_well_generation(0, n_points_potential, 0, 0.4);


    matched_from_line_density(GP->beta[0][0], GP->energy[0][0],
                              GP->charge, Beam->n_macroparticles,
                              Beam->dt, Beam->dE, GP->eta_0[0][0],
                              GP->t_rev[0], full_ring->fPotentialWell,
                              full_ring->fPotentialWellCoordinates,
                              line_density_opt);

}


#endif /* BEAMS_DISTRIBUTIONS_H_ */
