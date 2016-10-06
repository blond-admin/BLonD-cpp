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
#include <blond/utilities.h>
#include <blond/globals.h>
#include <blond/python.h>
#include <blond/trackers/Tracker.h>
#include <map>
#include <random>

void longitudinal_bigaussian(ftype sigma_dt, ftype sigma_dE = 0,
                             int seed = 0, bool reinsertion = false)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    if (GP->n_sections > 1) {
        std::cerr << "WARNING: longitudinal_bigaussian is not yet properly\n"
                  << "computed for several sections!\n";
    }
    if (RfP->n_rf > 1) {
        std::cerr << "WARNING: longitudinal_bigaussian for multiple RF\n"
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
    }

    // TODO if reinsertion == true
    if (reinsertion) {
        ;
    }
}


void matched_from_line_density(FullRingAndRf *full_ring,
                               std::map<std::string, std::string> line_density_opt,
                               std::string main_harmonic = "lowest_freq",
                               std::string plot = "",
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

    python::import();
    auto pFunc = python::import("distributions", "matched_from_line_density");
    auto pBeta = python::convert_double(GP->beta[0][0]);
    auto pEnergy = python::convert_double(GP->energy[0][0]);
    auto pCharge = python::convert_double(GP->charge);
    auto pNMacroparticles = python::convert_int(Beam->n_macroparticles);
    auto pEta0 = python::convert_double(GP->eta_0[0][0]);
    auto pTRev0 = python::convert_double(GP->t_rev[0]);
    auto pMainHarmonic = python::convert_string(main_harmonic);
    auto pPlot = python::convert_string(plot);
    auto pFigDir = python::convert_string(figdir);
    auto pHalfOption = python::convert_string(half_option);
    auto pNIterationsInput = python::convert_int(n_iterations_input);
    auto pSeed = python::convert_int(seed);
    auto pDT = python::convert_double_array(Beam->dt.data(), Beam->dt.size());
    auto pDE = python::convert_double_array(Beam->dE.data(), Beam->dE.size());

    auto pPotentialWell = python::convert_double_array(
                              full_ring->fPotentialWell.data(),
                              full_ring->fPotentialWell.size());

    auto pTimeCoord = python::convert_double_array(
                          full_ring->fPotentialWellCoordinates.data(),
                          full_ring->fPotentialWellCoordinates.size());

    auto pLineDensityOpt = python::convert_dictionary(line_density_opt);
    auto pExtraVoltageDict = python::convert_dictionary(extraVoltageDict);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pBeta, pEnergy,
                                            pCharge, pNMacroparticles,
                                            pDT, pDE, pEta0, pTRev0,
                                            pPotentialWell, pTimeCoord,
                                            pLineDensityOpt, pMainHarmonic,
                                            pPlot, pFigDir, pHalfOption,
                                            pExtraVoltageDict,
                                            pNIterationsInput, pSeed, NULL);

    if (!ret) {
        std::cerr << "[matched_from_line_density] An error occured while "
                  << "executing python code\n";
        exit(-1);
    }

}


void matched_from_distribution_density(FullRingAndRf *full_ring,
                                       std::map<std::string,
                                       std::string> distribution_opt,
                                       std::string main_harmonic = "lowest_freq",
                                       int n_iterations_input = 1,
                                       std::map<std::string, std::string> extraVoltageDict =
                                               std::map<std::string, std::string>(),
                                       int seed = 0
                                      )
{
    auto GP = Context::GP;
    auto Beam = Context::Beam;

    int n_points_potential = int(1e4);

    full_ring->potential_well_generation(0, n_points_potential, 0, 0.4);

    python::import();
    auto pFunc = python::import("distributions",
                                "matched_from_distribution_density");
    auto pBeta = python::convert_double(GP->beta[0][0]);
    auto pNMacroparticles = python::convert_int(Beam->n_macroparticles);
    auto pDT = python::convert_double_array(Beam->dt.data(), Beam->dt.size());
    auto pDE = python::convert_double_array(Beam->dE.data(), Beam->dE.size());
    auto pEnergy = python::convert_double(GP->energy[0][0]);
    auto pCharge = python::convert_double(GP->charge);
    auto pEta0 = python::convert_double(GP->eta_0[0][0]);
    auto pTRev0 = python::convert_double(GP->t_rev[0]);
    auto pMainHarmonic = python::convert_string(main_harmonic);
    auto pNIterationsInput = python::convert_int(n_iterations_input);
    auto pSeed = python::convert_int(seed);

    auto pPotentialWell = python::convert_double_array(
                              full_ring->fPotentialWell.data(),
                              full_ring->fPotentialWell.size());

    auto pTimeCoord = python::convert_double_array(
                          full_ring->fPotentialWellCoordinates.data(),
                          full_ring->fPotentialWellCoordinates.size());

    auto pDistributionOpt = python::convert_dictionary(distribution_opt);
    auto pExtraVoltageDict = python::convert_dictionary(extraVoltageDict);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pBeta, pEnergy,
                                            pCharge, pNMacroparticles,
                                            pDT, pDE, pEta0, pTRev0,
                                            pPotentialWell, pTimeCoord,
                                            pDistributionOpt, pMainHarmonic,
                                            pNIterationsInput,
                                            pExtraVoltageDict,
                                            pSeed, NULL);
    if (!ret) {
        std::cerr << "[matched_from_line_density] An error occured while "
                  << "executing python code\n";
        exit(-1);
    }

}


#endif /* BEAMS_DISTRIBUTIONS_H_ */
