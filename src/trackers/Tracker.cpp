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



RingAndRfSection::RingAndRfSection(RfParameters *rfp, solver_type solver,
                                   PhaseLoop *PhaseLoop, LHCNoiseFB *NoiseFB,
                                   bool periodicity, ftype dE_max,
                                   bool rf_kick_interp, Slices *Slices,
                                   TotalInducedVoltage *TotalInducedVoltage)
{
    auto GP = Context::GP;
    // auto Beam = Context::Beam;
    auto RfP = rfp;
    this->fRfP = RfP;

    this->elapsed_time = 0;
    this->solver = solver;
    this->PL = PhaseLoop;
    this->noiseFB = NoiseFB;
    this->periodicity = periodicity;
    this->dE_max = dE_max;
    this->rf_kick_interp = rf_kick_interp;
    this->slices = Slices;
    this->totalInducedVoltage = TotalInducedVoltage;

    this->acceleration_kick.resize(GP->n_turns + 1);
    for (uint i = 0; i < RfP->E_increment.size(); ++i)
        acceleration_kick[i] = -RfP->E_increment[i];

    if (solver != simple && solver != full) {
        std::cerr << "ERROR: Choice of longitudinal solver not recognized!\n"
                  << "Aborting...\n";
        exit(-1);
    }

    if (GP->alpha_order > 1) {
        solver = full;
    }

    if (rf_kick_interp && Slices == NULL) {
        std::cerr << "ERROR: A slices object is needed to use the"
                  << " kick_interp option\n";
    }

}


inline void RingAndRfSection::kick(const ftype *__restrict beam_dt,
                                   ftype *__restrict beam_dE,
                                   const int n_rf,
                                   const ftype *__restrict voltage,
                                   const ftype *__restrict omega_RF,
                                   const ftype *__restrict phi_RF,
                                   const int n_macroparticles,
                                   const ftype acc_kick)
{
    // KICK
    //#pragma omp parallel for collapse(2)
    for (int j = 0; j < n_rf; ++j) {
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; ++i) {
            const ftype a = omega_RF[j] * beam_dt[i] + phi_RF[j];
            beam_dE[i] += voltage[j] * mymath::fast_sin(a);
        }
    }

// SYNCHRONOUS ENERGY CHANGE
    #pragma omp parallel for
    for (int i = 0; i < n_macroparticles; ++i)
        beam_dE[i] += acc_kick;
}

// kick with periodicity

// inline void RingAndRfSection::kick(const ftype *__restrict beam_dt,
//                                    ftype *__restrict beam_dE,
//                                    const int n_rf,
//                                    const ftype *__restrict voltage,
//                                    const ftype *__restrict omega_RF,
//                                    const ftype *__restrict phi_RF,
//                                    const int n_macroparticles,
//                                    const ftype acc_kick,
//                                    const int_vector_t &filter)
// {
//     // KICK
//     //#pragma omp parallel for collapse(2)
//     for (int j = 0; j < n_rf; j++) {
//         for (const auto &i : filter) {
//             std::cout << i << "\n";
//             const ftype a = omega_RF[j] * beam_dt[i] + phi_RF[j];
//             beam_dE[i] += voltage[j] * mymath::fast_sin(a);
//         }
//     }

//     // SYNCHRONOUS ENERGY CHANGE
//     for (const auto &i : filter)
//         beam_dE[i] += acc_kick;
// }


inline void RingAndRfSection::drift(ftype *__restrict beam_dt,
                                    const ftype *__restrict beam_dE,
                                    const solver_type solver,
                                    const ftype T0,
                                    const ftype length_ratio,
                                    const uint alpha_order,
                                    const ftype eta_zero,
                                    const ftype eta_one,
                                    const ftype eta_two,
                                    const ftype beta,
                                    const ftype energy,
                                    const int n_macroparticles)
{

    const ftype T = T0 * length_ratio;

    if (solver == simple) {
        const ftype T_x_coeff = T * eta_zero / (beta * beta * energy);
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T_x_coeff * beam_dE[i];
    } else {
        const ftype coeff = 1. / (beta * beta * energy);
        const ftype eta0 = eta_zero * coeff;
        const ftype eta1 = eta_one * coeff * coeff;
        const ftype eta2 = eta_two * coeff * coeff * coeff;

        if (alpha_order == 1)
            #pragma omp parallel for
            for (int i = 0; i < n_macroparticles; i++)
                beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
        else if (alpha_order == 2)
            #pragma omp parallel for
            for (int i = 0; i < n_macroparticles; i++)
                beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i] -
                                         eta1 * beam_dE[i] * beam_dE[i]) -
                                   1.);
        else
            #pragma omp parallel for
            for (int i = 0; i < n_macroparticles; i++)
                beam_dt[i] +=
                    T * (1. / (1. - eta0 * beam_dE[i] -
                               eta1 * beam_dE[i] * beam_dE[i] -
                               eta2 * beam_dE[i] * beam_dE[i] * beam_dE[i]) -
                         1.);
    }
}

// drift with periodicity

// inline void RingAndRfSection::drift(ftype *__restrict beam_dt,
//                                     const ftype *__restrict beam_dE,
//                                     const solver_type solver,
//                                     const ftype T0,
//                                     const ftype length_ratio,
//                                     const uint alpha_order,
//                                     const ftype eta_zero,
//                                     const ftype eta_one,
//                                     const ftype eta_two,
//                                     const ftype beta,
//                                     const ftype energy,
//                                     const int n_macroparticles,
//                                     const int_vector_t &filter)
// {

//     ftype T = T0 * length_ratio;

//     if (solver == simple) {
//         ftype coeff = eta_zero / (beta * beta * energy);

//         for (const auto &i : filter)
//             beam_dt[i] += T * coeff * beam_dE[i];
//     }

//     else {
//         const ftype coeff = 1. / (beta * beta * energy);
//         const ftype eta0 = eta_zero * coeff;
//         const ftype eta1 = eta_one * coeff * coeff;
//         const ftype eta2 = eta_two * coeff * coeff * coeff;

//         if (alpha_order == 1)
//             for (const auto &i : filter)
//                 beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
//         else if (alpha_order == 2)
//             for (const auto &i : filter)
//                 beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i] -
//                                          eta1 * beam_dE[i] * beam_dE[i]) -
//                                    1.);
//         else
//             for (const auto &i : filter)
//                 beam_dt[i] +=
//                     T * (1. / (1. - eta0 * beam_dE[i] -
//                                eta1 * beam_dE[i] * beam_dE[i] -
//                                eta2 * beam_dE[i] * beam_dE[i] * beam_dE[i]) -
//                          1.);
//     }
// }

void RingAndRfSection::track()
{
    auto GP = Context::GP;
    auto RfP = fRfP;
    auto Beam = Context::Beam;

    if (!RfP->phi_noise.empty()) {
        if (noiseFB != NULL) {
            for (uint i = 0; i < RfP->phi_RF.size(); ++i)
                RfP->phi_RF[i][RfP->counter] +=
                    noiseFB->fX * RfP->phi_noise[i][RfP->counter];
        } else {
            for (uint i = 0; i < RfP->phi_RF.size(); ++i)
                RfP->phi_RF[i][RfP->counter] += RfP->phi_noise[i][RfP->counter];
        }
    }

    // Determine phase loop correction on RF phase and frequency
    if (PL != NULL && RfP->counter >= PL->delay)
        PL->track();

    if (periodicity) {
        // Change reference of all the particles on the right of the current
        // frame; these particles skip one kick and drift
        set_periodicity();
        const ftype tRev = GP->t_rev[RfP->counter + 1];

        if (!indices_right_outside.empty()) {
            f_vector_t insiders_dt, insiders_dE;
            insiders_dt.reserve(indices_inside_frame.size());
            insiders_dE.reserve(indices_inside_frame.size());
            for (const auto &i : indices_inside_frame) {
                insiders_dt.push_back(Beam->dt[i]);
                insiders_dE.push_back(Beam->dE[i]);
            }

            for (const auto &i : indices_right_outside)
                Beam->dt[i] -= tRev;

            // Synchronize the bunch with the particles that are on the right of
            // the current frame applying kick and drift to the bunch; after that
            // all the particle are in the new updated frame
            kick(insiders_dt, insiders_dE, RfP->counter);
            drift(insiders_dt, insiders_dE, RfP->counter + 1);
            int k = 0;
            for (const auto &i : indices_inside_frame) {
                Beam->dt[i] = insiders_dt[k];
                Beam->dE[i] = insiders_dE[k];
                k++;
            }

        } else {
            kick(Beam->dt, Beam->dE, RfP->counter);
            drift(Beam->dt, Beam->dE, RfP->counter + 1);
            // find left outside particles and kick, drift them one more time
            // indices_left_outside.clear();
            //#pragma omp parallel for reduction(+:a)
            for (uint i = 0; i < Beam->n_macroparticles; ++i) {
                if (Beam->dt[i] < 0)
                    indices_left_outside.push_back(i);
            }
        }

        if (!indices_left_outside.empty()) {
            f_vector_t left_dt, left_dE;
            left_dt.reserve(indices_left_outside.size());
            left_dE.reserve(indices_left_outside.size());
            for (const auto &i : indices_left_outside) {
                left_dt.push_back(Beam->dt[i] + tRev);
                left_dE.push_back(Beam->dE[i]);
            }

            kick(left_dt, left_dE, RfP->counter);
            drift(left_dt, left_dE, RfP->counter + 1);
            int k = 0;
            for (const auto &i : indices_left_outside) {
                Beam->dt[i] = left_dt[k];
                Beam->dE[i] = left_dE[k];
                k++;
            }
        }

        // std::cout << "right: " << indices_right_outside.size() << '\n';
        // std::cout << "inside: " << indices_inside_frame.size() << '\n';
        // std::cout << "left: " << indices_left_outside.size() << '\n';

    } else {
        if (rf_kick_interp) {
            // TODO complete & test this part
        } else {
            kick(Beam->dt, Beam->dE, RfP->counter);
        }
        drift(Beam->dt, Beam->dE, RfP->counter + 1);
    }

    if (dE_max > 0) horizontal_cut();

    RfP->counter++;
}


// TODO test this function
inline void RingAndRfSection::horizontal_cut()
{
    auto Beam = Context::Beam;

    uint i = 0;
    while (i < Beam->dE.size()) {
        if (Beam->dE[i] > -dE_max) {
            Beam->dE.erase(Beam->dE.begin() + i);
            Beam->dt.erase(Beam->dt.begin() + i);
            Beam->id.erase(Beam->id.begin() + i);
        } else {
            i++;
        }
    }
    Beam->n_macroparticles = Beam->dE.size();
}



RingAndRfSection::~RingAndRfSection() {}


void RingAndRfSection::set_periodicity()
{
    auto GP = Context::GP;
    auto RfP = fRfP;
    auto Beam = Context::Beam;

    indices_right_outside.clear();
    indices_inside_frame.clear();
    indices_left_outside.clear();
    // TODO I am not duplicating the insiders dE, dt
    // as done in the python version
    for (uint i = 0; i < Beam->n_macroparticles; ++i) {
        if (Beam->dt[i] > GP->t_rev[RfP->counter + 1])
            indices_right_outside.push_back(i);
        else if (Beam->dt[i] < GP->t_rev[RfP->counter + 1])
            indices_inside_frame.push_back(i);
    }
}

void RingAndRfSection::kick(f_vector_t &beam_dt, f_vector_t &beam_dE, const uint index)
{
    auto RfP = fRfP;
    // auto Beam = Context::Beam;

    auto vol = new ftype[RfP->n_rf];
    auto omeg = new ftype[RfP->n_rf];
    auto phi = new ftype[RfP->n_rf];

    for (uint i = 0; i < RfP->n_rf; ++i) {
        vol[i] = RfP->voltage[i][index];
        omeg[i] = RfP->omega_RF[i][index];
        phi[i] = RfP->phi_RF[i][index];
    }

    kick(beam_dt.data(), beam_dE.data(), RfP->n_rf, vol, omeg, phi,
         beam_dt.size(), acceleration_kick[index]);

    delete[] vol;
    delete[] omeg;
    delete[] phi;
}

// void RingAndRfSection::kick(const int_vector_t &filter, const uint index)
// {
//     auto RfP = fRfP;
//     auto Beam = Context::Beam;

//     auto vol = new ftype[RfP->n_rf];
//     auto omeg = new ftype[RfP->n_rf];
//     auto phi = new ftype[RfP->n_rf];

//     for (uint i = 0; i < RfP->n_rf; ++i) {
//         vol[i] = RfP->voltage[i][index];
//         omeg[i] = RfP->omega_RF[i][index];
//         phi[i] = RfP->phi_RF[i][index];
//     }

//     kick(Beam->dt.data(), Beam->dE.data(), RfP->n_rf, vol, omeg, phi,
//          Beam->n_macroparticles, acceleration_kick[index], filter);

//     delete[] vol;
//     delete[] omeg;
//     delete[] phi;
// }

// void RingAndRfSection::drift(const int_vector_t &filter, const uint index)
// {
//     auto GP = Context::GP;
//     auto RfP = fRfP;
//     auto Beam = Context::Beam;

//     drift(Beam->dt.data(), Beam->dE.data(), solver, GP->t_rev[index],
//           RfP->length_ratio, GP->alpha_order, RfP->eta_0(index),
//           RfP->eta_1(index), RfP->eta_2(index), RfP->beta(index),
//           RfP->energy(index), Beam->n_macroparticles, filter);
// }

void RingAndRfSection::drift(f_vector_t &beam_dt, f_vector_t &beam_dE, const uint index)
{
    auto GP = Context::GP;
    auto RfP = fRfP;
    // auto Beam = Context::Beam;

    drift(beam_dt.data(), beam_dE.data(), solver, GP->t_rev[index],
          RfP->length_ratio, GP->alpha_order, RfP->eta_0(index),
          RfP->eta_1(index), RfP->eta_2(index), RfP->beta(index),
          RfP->energy(index), beam_dt.size());
}

FullRingAndRf::FullRingAndRf(std::vector<RingAndRfSection *> &RingList)
{
    fRingList = RingList;

    fRingCircumference = 0;
    for (auto &ring : fRingList)
        fRingCircumference += ring->fRfP->section_length;

    fRingRadius = fRingCircumference / (2 * constant::pi);
}

FullRingAndRf::~FullRingAndRf() {}

void FullRingAndRf::track()
{
    // Loops over all the RingAndRFSection.track methods.
    for (auto &ring : fRingList)
        ring->track();
}

void FullRingAndRf::potential_well_generation(const uint turn,
        const uint n_points,
        const ftype option,
        const ftype dt_margin_percent)
{
    auto GP = Context::GP;
    f_vector_t voltages;
    f_vector_t omega_rf;
    f_vector_t phi_offsets;

    for (const auto &ring : fRingList) {
        for (uint i = 0; i < ring->fRfP->n_rf; ++i) {
            voltages.push_back(ring->fRfP->voltage[i][turn]);
            omega_rf.push_back(ring->fRfP->omega_RF[i][turn]);
            phi_offsets.push_back(ring->fRfP->phi_RF[i][turn]);
        }
    }

    ftype main_omega_rf;
    if (option == 0) {
        auto k = mymath::min(omega_rf.data(), omega_rf.size());
        main_omega_rf = omega_rf[k];
    } else if (option == 1) {
        auto k = mymath::max(voltages.data(), voltages.size());
        auto maxV = voltages[k];

        f_vector_t temp;

        for (uint i = 0; i < voltages.size(); ++i)
            if (voltages[i] == maxV)
                temp.push_back(omega_rf[i]);

        k = mymath::min(temp.data(), temp.size());
        main_omega_rf = temp[k];
    } else {
        f_vector_t temp;
        std::copy_if(omega_rf.begin(), omega_rf.end(), back_inserter(temp),
        [option](const ftype x) { return x == option; });
        if (temp.empty()) {
            std::cerr << "[ERROR] The desired harmonic to compute"
                      << "the potential well does not"
                      << "match the RF parameters...\n";
            exit(-1);
        }
        main_omega_rf = option;//temp[k];
    }

    auto time_array_margin =
        dt_margin_percent * 2 * constant::pi / main_omega_rf;

    auto slippage_factor = fRingList[0]->fRfP->eta_0(turn);

    auto first_dt = -time_array_margin / 2;
    auto last_dt = 2 * constant::pi / main_omega_rf + time_array_margin / 2;

    f_vector_t time_array(n_points);
    mymath::linspace(time_array.data(), first_dt, last_dt, n_points);

    fTotalVoltage.resize(time_array.size());

    for (int i = 0; i < (int)time_array.size(); ++i) {
        ftype sum = 0.0;
        for (int j = 0; j < (int)voltages.size(); ++j) {
            sum += voltages[j] * mymath::fast_sin(omega_rf[j] * time_array[i] +
                                                  phi_offsets[j]);
        }
        fTotalVoltage[i] = sum;
    }

    const auto eom_factor_potential =
        mymath::sign(slippage_factor) * GP->charge / GP->t_rev[turn];

    auto tempVec = fTotalVoltage;
    for (int i = 0; i < (int)fTotalVoltage.size(); ++i)
        tempVec[i] =
            eom_factor_potential *
            (tempVec[i] -
             (-fRingList[0]->acceleration_kick[turn]) / std::abs(GP->charge));

    auto tempArr = mymath::cum_trapezoid(
                       tempVec.data(), time_array[1] - time_array[0], n_points);

    std::transform(tempArr.begin(), tempArr.end(), tempArr.begin(),
    [](ftype x) { return -x; });

    auto min_i = mymath::min(tempArr.data(), tempArr.size());
    auto min = tempArr[min_i] > 0 ? 0 : tempArr[min_i];

    fPotentialWell.resize(tempArr.size() + 1);
    fPotentialWell[0] = - min;

    for (uint i = 0; i < tempArr.size(); ++i)
        fPotentialWell[i + 1] = tempArr[i] - min;

    fPotentialWellCoordinates = time_array;
}