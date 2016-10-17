/*
 * RfParameters.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */


#include <blond/input_parameters/RfParameters.h>

/*
 :How to use RF programs:

 - For 1 RF system and constant values of V, h or phi, just input the single
 value
 - For 1 RF system and varying values of V, h or phi, input an array of n_turns
 values
 - For several RF systems and constant values of V, h or phi, input lists of
 single values
 - For several RF systems and varying values of V, h or phi, input lists of
 arrays of n_turns values
 */


double RfParameters::eta_tracking(const Beams *beam, const int counter,
                                  const double dE) const
{

    double eta = 0;
    if (alpha_order == 1)
        eta = eta_0[counter];
    else {
        double delta =
            dE / ((beam->beta) * (beam->beta) * beam->energy);
        eta += eta_0[counter] * 1;
        if (alpha_order > 0)
            eta += eta_1[counter] * delta;
        if (alpha_order > 1)
            eta += eta_2[counter] * delta * delta;
        if (alpha_order > 2)
            std::cerr << "WARNING: Momentum compaction factor is implemented"
                      << "only up to 2nd order\n";
    }
    return eta;
}



f_vector_t RfParameters::calc_phi_s(RfParameters *rfp,
                                    acc_sys_t acc_sys)
{
    /*
     | *The synchronous phase calculated from the rate of momentum change.*
     | *Below transition, for decelerating bucket: phi_s is in (-Pi/2,0)*
     | *Below transition, for accelerating bucket: phi_s is in (0,Pi/2)*
     | *Above transition, for accelerating bucket: phi_s is in (Pi/2,Pi)*
     | *Above transition, for decelerating bucket: phi_s is in (Pi,3Pi/2)*
     | *The synchronous phase is calculated at a certain moment.*
     | *Uses beta, energy averaged over the turn.*
     */

    auto n_turns = rfp->n_turns;
    auto n_rf = rfp->n_rf;
    // std::cout << "n_turns: " << n_turns << "\n";
    f_vector_t out(n_turns + 1);
    // double eta0 = rf_params->eta0;
    if (acc_sys == RfParameters::acc_sys_t::as_single) {

        f_vector_t denergy = rfp->E_increment;
        denergy.push_back(rfp->E_increment.back());

        f_vector_t acceleration_ratio(n_turns + 1);
        for (int i = 0; i < n_turns + 1; ++i)
            acceleration_ratio[i] =
                denergy[i] / (rfp->charge * rfp->voltage[rfp->section_index][i]);

        for (int i = 0; i < n_turns + 1; ++i)
            if (acceleration_ratio[i] > 1 || acceleration_ratio[i] < -1)
                dprintf("Warning!!! Acceleration is not possible (momentum "
                        "increment "
                        "is too big or voltage too low) at index %d\n",
                        i);

        for (int i = 0; i < n_turns + 1; ++i)
            out[i] = asin(acceleration_ratio[i]);

        // double *eta0_middle_points = new double[n_turns +1];
        for (int i = 0; i < n_turns; ++i) {
            double middle = (rfp->eta_0[i] + rfp->eta_0[i + 1]) / 2;
            if (middle > 0)
                out[i] = constant::pi - out[i];
            else
                out[i] = constant::pi + out[i];
        }
        if (rfp->eta_0[n_turns] > 0)
            out[n_turns] = constant::pi - out[n_turns];
        else
            out[n_turns] = constant::pi + out[n_turns];

        return out;

    } else if (acc_sys == RfParameters::acc_sys_t::all) {
        /*
         In this case, all the rf systems are accelerating, phi_s is
         calculated accordingly with respect to the fundamental frequency (the
         minimum
         of the potential well is taken)
         */

        f_vector_t transition_phase_offset(n_turns + 1);

        for (int i = 0; i < n_turns + 1; ++i) {
            out[i] = 0;
            if (rfp->eta_0[i] > 0)
                transition_phase_offset[i] = constant::pi;
            else
                transition_phase_offset[i] = 0;
        }
        f_vector_t phase_array(1000);
        mymath::linspace(phase_array.data(), -constant::pi * 1.2,
                         constant::pi * 1.2, 1000);

        for (int i = 0; i < n_turns; ++i) {
            f_vector_t totalrf(1000, 0); // = { };
            for (int j = 0; j < n_rf; ++j) {
                double min = rfp->harmonic[0][i];
                for (int k = 0; k < n_rf; ++k)
                    min = std::min(min, rfp->harmonic[k][i]);
                for (int k = 0; k < 1000; ++k) {
                    totalrf[k] +=
                        rfp->voltage[j][i + 1] *
                        std::sin((rfp->harmonic[j][i + 1] / min) *
                                 (phase_array[k] +
                                  transition_phase_offset[i + 1]) +
                                 rfp->phi_offset[j][i + 1]);
                }
            }
            int transition_factor = transition_phase_offset[i] == 0 ? +1 : -1;

            double potential_well[1000] = {0};
            f_vector_t f(1000);
            for (int k = 0; k < 1000; ++k) {
                f[k] = totalrf[k] - rfp->E_increment[i] / rfp->charge;
            }

            auto trap = mymath::cum_trapezoid(f.data(),
                                              phase_array[1] - phase_array[0],
                                              1000);

            for (int k = 0; k < 1000; ++k) {
                potential_well[k] = transition_factor * trap[k];
            }

            // TODO is this correct? line BLonD-minimal::rf_parameter.py:334
            out[i + 1] = phase_array[mymath::min(potential_well, 1000, 1)] +
                         transition_phase_offset[i + 1];
        }
        out[0] = out[1];
        return out;

    } else if (acc_sys == RfParameters::acc_sys_t::first) {
        /*
         Only the first rf system is accelerating, so we have to correct the
         phi_offset of the other rf_systems such that p_increment relates
         only to the first rf
         */
        ;
    } else {
        std::cerr << "Did not recognize the option accelerating_systems in calc_phi_s "
                  << "function\n";
        exit(-1);
    }

    // TODO why n_turns and not n_turns+1?
    if (rfp->eta_0[0] > 0.)
        out.resize(n_turns, constant::pi);
    else
        out.resize(n_turns, 0.);
    return out;
}
