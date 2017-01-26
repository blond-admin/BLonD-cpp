/*
 * Tracker.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include <blond/constants.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <iterator>
#include <blond/vector_math.h>

using namespace blond;
using namespace std;
using namespace mymath;


inline void RingAndRfSection::kick(const double *__restrict beam_dt,
                                   double *__restrict beam_dE,
                                   const int n_rf,
                                   const double *__restrict voltage,
                                   const double *__restrict omega_rf,
                                   const double *__restrict phi_rf,
                                   const int n_macroparticles,
                                   const double acc_kick)
{
    // KICK
    //#pragma omp parallel for collapse(2)
    for (int j = 0; j < n_rf; ++j) {
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; ++i) {
            const double a = omega_rf[j] * beam_dt[i] + phi_rf[j];
            beam_dE[i] += voltage[j] * fast_sin(a);
        }
    }

// SYNCHRONOUS ENERGY CHANGE
    #pragma omp parallel for
    for (int i = 0; i < n_macroparticles; ++i)
        beam_dE[i] += acc_kick;
}


inline void RingAndRfSection::drift(double *__restrict beam_dt,
                                    const double *__restrict beam_dE,
                                    const solver_type solver,
                                    const double T0,
                                    const double length_ratio,
                                    const int alpha_order,
                                    const double eta_zero,
                                    const double eta_one,
                                    const double eta_two,
                                    const double beta,
                                    const double energy,
                                    const int n_macroparticles)
{

    const double T = T0 * length_ratio;

    if (solver == simple) {
        const double T_x_coeff = T * eta_zero / (beta * beta * energy);
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T_x_coeff * beam_dE[i];
    } else {
        const double coeff = 1. / (beta * beta * energy);
        const double eta0 = eta_zero * coeff;
        const double eta1 = eta_one * coeff * coeff;
        const double eta2 = eta_two * coeff * coeff * coeff;
        if (alpha_order == 1) {
            #pragma omp parallel for
            for (int i = 0; i < n_macroparticles; i++)
                beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
        } else if (alpha_order == 2) {
            #pragma omp parallel for
            for (int i = 0; i < n_macroparticles; i++)
                beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i] -
                                         eta1 * beam_dE[i] * beam_dE[i]) -
                                   1.);
        } else {
            #pragma omp parallel for
            for (int i = 0; i < n_macroparticles; i++)
                beam_dt[i] +=
                    T * (1. / (1. - eta0 * beam_dE[i] -
                               eta1 * beam_dE[i] * beam_dE[i] -
                               eta2 * beam_dE[i] * beam_dE[i] * beam_dE[i]) -
                         1.);
        }
    }
}



void RingAndRfSection::track()
{

    if (!phi_noise.empty()) {
        if (noiseFB != NULL) {
            for (uint i = 0; i < phi_rf.size(); ++i)
                phi_rf[i][counter] +=
                    noiseFB->fX * phi_noise[i][counter];
        } else {
            for (uint i = 0; i < phi_rf.size(); ++i)
                phi_rf[i][counter] += phi_noise[i][counter];
        }
    }

    // Determine phase loop correction on RF phase and frequency
    if (PL != NULL && counter >= (int) PL->delay)
        PL->track();

    if (periodicity) {
        // Change reference of all the particles on the right of the current
        // frame; these particles skip one kick and drift
        set_periodicity();
        const double tRev = t_rev[counter + 1];

        if (!indices_right_outside.empty()) {
            f_vector_t insiders_dt, insiders_dE;
            insiders_dt.reserve(indices_inside_frame.size());
            insiders_dE.reserve(indices_inside_frame.size());
            for (const auto &i : indices_inside_frame) {
                insiders_dt.push_back(beam->dt[i]);
                insiders_dE.push_back(beam->dE[i]);
            }

            for (const auto &i : indices_right_outside)
                beam->dt[i] -= tRev;

            // Synchronize the bunch with the particles that are on the right of
            // the current frame applying kick and drift to the bunch; after that
            // all the particle are in the new updated frame
            kick(insiders_dt, insiders_dE, counter);
            drift(insiders_dt, insiders_dE, counter + 1);
            int k = 0;
            for (const auto &i : indices_inside_frame) {
                beam->dt[i] = insiders_dt[k];
                beam->dE[i] = insiders_dE[k];
                k++;
            }

        } else {
            kick(beam->dt, beam->dE, counter);
            drift(beam->dt, beam->dE, counter + 1);
            // find left outside particles and kick, drift them one more time
            // indices_left_outside.clear();
            //#pragma omp parallel for reduction(+:a)
            for (int i = 0; i < beam->n_macroparticles; ++i) {
                if (beam->dt[i] < 0)
                    indices_left_outside.push_back(i);
            }
        }

        if (!indices_left_outside.empty()) {
            f_vector_t left_dt, left_dE;
            left_dt.reserve(indices_left_outside.size());
            left_dE.reserve(indices_left_outside.size());
            for (const auto &i : indices_left_outside) {
                left_dt.push_back(beam->dt[i] + tRev);
                left_dE.push_back(beam->dE[i]);
            }

            kick(left_dt, left_dE, counter);
            drift(left_dt, left_dE, counter + 1);
            int k = 0;
            for (const auto &i : indices_left_outside) {
                beam->dt[i] = left_dt[k];
                beam->dE[i] = left_dE[k];
                k++;
            }
        }

        // cout << "right: " << indices_right_outside.size() << '\n';
        // cout << "inside: " << indices_inside_frame.size() << '\n';
        // cout << "left: " << indices_left_outside.size() << '\n';

    } else {
        if (rf_kick_interp) {
            // TODO test this part
            rf_voltage_calculation(counter, slices);
            fTotalVoltage = fRfVoltage;

            if (totalInducedVoltage != NULL)
                fTotalVoltage += totalInducedVoltage->fInducedVoltage;

            fRfVoltage *= charge;

            linear_interp_kick(beam->dt.data(), beam->dE.data(),
                               fRfVoltage.data(), slices->bin_centers.data(),
                               slices->n_slices, beam->n_macroparticles);
        } else {
            kick(beam->dt, beam->dE, counter);
        }
        drift(beam->dt, beam->dE, counter + 1);
    }

    if (dE_max > 0) horizontal_cut();

    counter++;
}


// TODO test this function
inline void RingAndRfSection::horizontal_cut()
{
    uint i = 0;
    while (i < beam->dE.size()) {
        if (beam->dE[i] > -dE_max) {
            beam->dE.erase(beam->dE.begin() + i);
            beam->dt.erase(beam->dt.begin() + i);
            beam->id.erase(beam->id.begin() + i);
        } else {
            i++;
        }
    }
    beam->n_macroparticles = beam->dE.size();
}

void RingAndRfSection::rf_voltage_calculation(int turn, Slices *slices)
{
    // Calculating the RF voltage seen by the beam at a given turn,
    // needs a Slices object.

    auto vol = new double[n_rf];
    auto omeg = new double[n_rf];
    auto phi = new double[n_rf];

    for (int i = 0; i < n_rf; ++i) {
        vol[i] = voltage[i][turn];
        omeg[i] = omega_rf[i][turn];
        phi[i] = phi_rf[i][turn];
    }

    fRfVoltage.resize(slices->bin_centers.size());

    for (uint j = 0; j < slices->bin_centers.size(); j++) {
        double sum = 0.0;
        for (int i = 0; i < n_rf; i++) {
            sum += vol[i] * sin(omeg[i] * slices->bin_centers[j] + phi[i]);
        }
        fRfVoltage[j] = sum;
    }


    delete[] vol;
    delete[] omeg;
    delete[] phi;
}


void RingAndRfSection::set_periodicity()
{

    indices_right_outside.clear();
    indices_inside_frame.clear();
    indices_left_outside.clear();
    // TODO I am not duplicating the insiders dE, dt
    // as done in the python version
    for (int i = 0; i < beam->n_macroparticles; ++i) {
        if (beam->dt[i] > t_rev[counter + 1])
            indices_right_outside.push_back(i);
        else if (beam->dt[i] < t_rev[counter + 1])
            indices_inside_frame.push_back(i);
    }
}

void RingAndRfSection::kick(f_vector_t &beam_dt, f_vector_t &beam_dE,
                            const int index)
{

    auto vol = new double[n_rf];
    auto omeg = new double[n_rf];
    auto phi = new double[n_rf];

    for (int i = 0; i < n_rf; ++i) {
        vol[i] = voltage[i][index];
        omeg[i] = omega_rf[i][index];
        phi[i] = phi_rf[i][index];
    }

    kick(beam_dt.data(), beam_dE.data(), n_rf, vol, omeg, phi,
         beam_dt.size(), acceleration_kick[index]);

    delete[] vol;
    delete[] omeg;
    delete[] phi;
}


void RingAndRfSection::drift(f_vector_t &beam_dt, f_vector_t &beam_dE,
                             const int index)
{

    drift(beam_dt.data(), beam_dE.data(), solver, t_rev[index],
          length_ratio, alpha_order, eta_0[index],
          eta_1[index], eta_2[index], rfp->beta[index],
          rfp->energy[index], beam_dt.size());
}

FullRingAndRf::FullRingAndRf(const vector<RingAndRfSection *> &RingList)
{
    fRingList = RingList;

    fRingCircumference = 0;
    for (auto &ring : fRingList)
        fRingCircumference += ring->section_length;

    fRingRadius = fRingCircumference / (2 * constant::pi);
}

FullRingAndRf::~FullRingAndRf() {}

void FullRingAndRf::track()
{
    // Loops over all the RingAndRFSection.track methods.
    for (auto &ring : fRingList)
        ring->track();
}

void FullRingAndRf::potential_well_generation(const int turn,
        const int n_points,
        const double option,
        const double dt_margin_percent)
{
    f_vector_t voltages;
    f_vector_t omega_rf;
    f_vector_t phi_offsets;
    auto charge = fRingList.back()->charge;
    auto &t_rev = fRingList.back()->t_rev;
    for (const auto &ring : fRingList) {
        for (int i = 0; i < ring->n_rf; ++i) {
            voltages.push_back(ring->voltage[i][turn]);
            omega_rf.push_back(ring->omega_rf[i][turn]);
            phi_offsets.push_back(ring->phi_rf[i][turn]);
        }
    }

    double main_omega_rf;
    if (option == 0) {
        main_omega_rf = *min_element(ALL(omega_rf));
    } else if (option == 1) {
        auto maxV = *min_element(ALL(voltages));
        f_vector_t temp;

        for (uint i = 0; i < voltages.size(); ++i)
            if (voltages[i] == maxV)
                temp.push_back(omega_rf[i]);

        main_omega_rf = *min_element(ALL(omega_rf));
    } else {
        f_vector_t temp;
        copy_if(omega_rf.begin(), omega_rf.end(), back_inserter(temp),
        [option](const double x) { return x == option; });
        if (temp.empty()) {
            cerr << "[ERROR] The desired harmonic to compute"
                 << "the potential well does not"
                 << "match the RF parameters...\n";
            exit(-1);
        }
        main_omega_rf = option;//temp[k];
    }

    // cout << "main_omega_rf: " << main_omega_rf << "\n";

    double time_array_margin =
        dt_margin_percent * 2 * constant::pi / main_omega_rf;

    double slippage_factor = fRingList[0]->eta_0[turn];

    double first_dt = -time_array_margin / 2;
    double last_dt = 2 * constant::pi / main_omega_rf + time_array_margin / 2;

    f_vector_t time_array = linspace(first_dt, last_dt, n_points);

    fTotalVoltage.resize(time_array.size());

    for (int i = 0; i < (int)time_array.size(); ++i) {
        double sum = 0.0;
        for (int j = 0; j < (int)voltages.size(); ++j) {
            sum += voltages[j] * fast_sin(omega_rf[j] * time_array[i] +
                                          phi_offsets[j]);
        }
        fTotalVoltage[i] = sum;
    }

    const double eom_factor_potential = sign(slippage_factor) * charge
                                        / t_rev[turn];

    // cout << "eom_factor_potential: " << eom_factor_potential << "\n";
    // cout << "fTotalVoltage: " << fTotalVoltage << "\n";
    // cout << "time_array: " << time_array[1] - time_array[0] << "\n";

    fPotentialWell = 0.0 - cum_trapezoid(eom_factor_potential / abs(charge)
                                         * (fTotalVoltage
                                            + fRingList[0]->acceleration_kick[turn]),
                                         time_array[1] - time_array[0]);


    fPotentialWell.insert(fPotentialWell.begin(), 0);
    // cout << "fPotentialWell: " << fPotentialWell;
    fPotentialWell -= *min_element(ALL(fPotentialWell));

    fPotentialWellCoordinates = time_array;
}