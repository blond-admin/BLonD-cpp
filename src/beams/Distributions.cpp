/*
 * Beams.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */
#include <blond/beams/Distributions.h>
#include <string>
#include <blond/python.h>
#include <iostream>

using namespace std;

// void matched_from_line_density(Beams *beam,
//                                FullRingAndRf *full_ring,
//                                map<string, string> line_density_opt,
//                                FullRingAndRf::main_harmonic_t main_harmonic_option,
//                                TotalInducedVoltage *totVolt,
//                                string plot,
//                                string figdir,
//                                string half_option,
//                                map<string, f_vector_t> extraVoltageDict,
//                                int n_iterations_input,
//                                int seed
//                               )
// {
//     // NOTE seed random engine
//     // not setting line_density_opt["exponent"] to null
//     double slippage_factor = full_ring->fRingList[0]->eta_0[0];
//     double eom_factor_dE = abs(slippage_factor) / (2 * beam->beta
//                            * beam->beta * beam->energy);

//     double eom_factor_potential = mymath::sign(slippage_factor)
//                                   * beam->charge
//                                   / full_ring->fRingList[0]->t_rev[0];
//     int n_points_potential = 1e4;
//     full_ring->potential_well_generation(0, n_points_potential,
//                                          main_harmonic_option,
//                                          0.4);

//     auto potential_well_array = full_ring->fPotentialWell;
//     auto time_coord_array = full_ring->fPotentialWellCoordinates;

//     f_vector_t extra_potential;
//     int n_points_line_den = 0;
//     f_vector_t line_density, time_line_den;
//     // NOTE what happens when extraVoltageDict is not empty??
//     double line_den_resolution = 0.;

//     if (!extraVoltageDict.empty()) {
//         auto &extra_voltage_time_input = extraVoltageDict["time_array"];
//         auto &extra_voltage_input = extraVoltageDict["voltage_array"];
//         auto extra_potential_input = mymath::cum_trapezoid(
//                                          extra_voltage_input.data(),
//                                          extra_voltage_input[1]
//                                          - extra_voltage_input[0],
//                                          extra_voltage_input.size()
//                                      );
//         extra_potential_input.insert(extra_potential_input.begin(), 0);
//         for (auto &e : extra_potential_input)
//             e *= - eom_factor_potential;
//         mymath::lin_interp(time_coord_array, extra_voltage_time_input,
//                            extra_potential_input, extra_potential);
//     }
//     if (line_density_opt.find("type") == line_density_opt.end()) {
//         cerr << "[matched_from_line_density] The input for the"
//              << "matched_from_line_density function was not recognized\n";
//         exit(-1);
//     }

//     if (line_density_opt["type"] != "user_input") {
//         n_points_line_den = 1e4;
//         time_line_den = mymath::linspace(time_coord_array.front(),
//                                          time_coord_array.back(),
//                                          n_points_line_den);
//         line_den_resolution = time_line_den[1] - time_line_den[0];
//         auto exponent = (line_density_opt.find("exponent") == line_density_opt.end()) ?
//                         0. : stod(line_density_opt["exponent"]);
//         line_density = line_density_function(time_line_den,
//                                              line_density_opt["type"],
//                                              stod(line_density_opt["bunch_length"]),
//                                              (time_coord_array.front()
//                                               - time_coord_array.back()) / 2.,
//                                              exponent);

//         auto min = *min_element(line_density.begin(), line_density.end());
//         for (auto &l : line_density) l -= min;
//         const auto sum = accumulate(line_density.begin(), line_density.end(), 0.);
//         for (auto &l : line_density) l = l / sum * beam->n_macroparticles;

//     } else { // (line_density_opt["type"] == "user_input") {
//         time_line_den = util::string_to_double_vector(line_density_opt["time_line_den"]);
//         n_points_line_den = time_line_den.size();
//         line_den_resolution = time_line_den[1] - time_line_den[0];

//         line_density = util::string_to_double_vector(
//                            line_density_opt["line_density"]);

//         auto min = *min_element(line_density.begin(), line_density.end());
//         for (auto &l : line_density) l -= min;
//         const auto sum = accumulate(line_density.begin(), line_density.end(), 0.);
//         for (auto &l : line_density) l = l / sum * beam->n_macroparticles;

//     }

//     f_vector_t induced_potential_final, induced_potential, time_induced_voltage;
//     int n_iterations = 1;

//     if (totVolt != nullptr) {
//         // NOTE continue here
//         auto rfp = totVolt->fSlices->rfp;
//         auto slices = Slices(rfp, beam, n_points_line_den);
//         slices.n_macroparticles = line_density;
//         slices.bin_centers = time_line_den;
//         slices.edges = mymath::linspace(time_line_den[0]
//                                         - (time_line_den[1] - time_line_den[0]) / 2,
//                                         time_line_den.back()
//                                         + (time_line_den[1] - time_line_den[0]) / 2,
//                                         time_line_den.size() + 1);
//         // fit option is already normal
//         auto induced_voltage_object = *totVolt;
//         // cout << "totVolt slices[0]: " << totVolt->fSlices->bin_centers[0] << "\n";
//         // cout << "object slices[0]: " << induced_voltage_object.fSlices->bin_centers[0] << "\n";
//         induced_voltage_object.reprocess(&slices);
//         // cout << "totVolt slices[0]: " << totVolt->fSlices->bin_centers[0] << "\n";
//         // cout << "object slices[0]: " << induced_voltage_object.fSlices->bin_centers[0] << "\n";

//         int induced_voltage_length = 1.5 * n_points_line_den;
//         auto induced_voltage = induced_voltage_object.induced_voltage_sum(
//                                    beam, induced_voltage_length);
//         time_induced_voltage = mymath::linspace(time_line_den[0],
//                                                 time_line_den[0] + (induced_voltage_length - 1)
//                                                 * line_den_resolution,
//                                                 induced_voltage_length);

//         induced_potential = mymath::cum_trapezoid(induced_voltage,
//                             time_induced_voltage[1] - time_induced_voltage[0]);
//         induced_potential.insert(induced_voltage.begin(), 0);
//         n_iterations = n_iterations_input;

//     }

//     for (int i = 0; i < n_iterations; i++) {
//         if (totVolt != nullptr)
//             mymath::lin_interp(time_coord_array, time_induced_voltage,
//                                induced_potential, induced_potential_final);
//         // TODO is there a need for all that or size of potential_well_array,
//         // induced_potential_final and total_potential is always the same?
//         uint size = max(potential_well_array.size(),
//                         max(induced_potential_final.size(),
//                             extra_potential.size()));

//         f_vector_t total_potential(size, 0.);

//         for (uint j = 0; j < size; j++) {
//             if (j < potential_well_array.size())
//                 total_potential[j] += potential_well_array[j];
//             if (j < induced_potential_final.size())
//                 total_potential[j] += induced_potential_final[j];
//             if (j < extra_potential.size())
//                 total_potential[j] += extra_potential[j];
//         }

//         //TODO I need minmax_location and potential_well_cut :(


//     }

// }


void matched_from_line_density(Beams *beam,
                               FullRingAndRf *full_ring,
                               map<string, string> line_density_opt,
                               string main_harmonic,
                               TotalInducedVoltage *totVolt,
                               string plot,
                               string figdir,
                               string half_option,
                               map<string, string> extraVoltageDict,
                               int n_iterations_input,
                               int seed
                              )
{
    auto GP = Context::GP;

    int n_points_potential = int(1e4);

    full_ring->potential_well_generation(0, n_points_potential, 0, 0.4);

    python::import();
    auto pFunc = python::import("distributions", "matched_from_line_density");
    auto pBeta = python::convert_double(beam->beta);
    auto pEnergy = python::convert_double(beam->energy);
    auto pCharge = python::convert_double(beam->charge);
    auto pNMacroparticles = python::convert_int(beam->n_macroparticles);
    auto pEta0 = python::convert_double(GP->eta_0[0][0]);
    auto pTRev0 = python::convert_double(GP->t_rev[0]);
    auto pMainHarmonic = python::convert_string(main_harmonic);
    auto pPlot = python::convert_string(plot);
    auto pFigDir = python::convert_string(figdir);
    auto pHalfOption = python::convert_string(half_option);
    auto pNIterationsInput = python::convert_int(n_iterations_input);
    auto pSeed = python::convert_int(seed);
    auto pDT = python::convert_double_array(beam->dt.data(), beam->dt.size());
    auto pDE = python::convert_double_array(beam->dE.data(), beam->dE.size());

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
                                       map<string, string> distribution_opt,
                                       string main_harmonic,
                                       int n_iterations_input,
                                       map<string, string> extraVoltageDict,
                                       int seed
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

    int counter = RfP->counter;
    ftype harmonic = RfP->harmonic[0][counter];
    ftype energy = GP->energy[0][counter];
    ftype beta = GP->beta[0][counter];
    ftype omega_rf = RfP->omega_rf[0][counter];
    ftype phi_s = RfP->phi_s[counter];
    ftype phi_rf = RfP->phi_rf[0][counter];
    ftype eta0 = RfP->eta_0[counter];

    if (sigma_dE == 0) {
        auto voltage = GP->charge * RfP->voltage[0][counter];
        auto phi_b = omega_rf * sigma_dt + phi_s;
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
        for (int i = 0; i < Beam->n_macroparticles; ++i) {
            ftype r = random[i % random.size()];
            if (eta0 > 0)
                Beam->dt[i] = sigma_dt * r
                              + (phi_s - phi_rf) / omega_rf;
            else
                Beam->dt[i] = sigma_dt * r
                              + (phi_s - phi_rf - constant::pi) / omega_rf;
            Beam->dE[i] = sigma_dE * r;
        }
    } else {
        default_random_engine generator(seed);
        normal_distribution<ftype> distribution(0.0, 1.0);
        for (int i = 0; i < Beam->n_macroparticles; ++i) {
            ftype r = distribution(generator);
            if (eta0 > 0)
                Beam->dt[i] = sigma_dt * r + (phi_s - phi_rf) / omega_rf;
            else
                Beam->dt[i] =
                    sigma_dt * r + (phi_s - phi_rf - constant::pi) / omega_rf;
            r = distribution(generator);
            Beam->dE[i] = sigma_dE * r;
        }
    }

    // TODO if reinsertion == true
    if (reinsertion) {
        ;
    }
}