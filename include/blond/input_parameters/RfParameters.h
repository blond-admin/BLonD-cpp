/*
 * RfParameters.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef INPUT_PARAMETERS_RFPARAMETERS_H_
#define INPUT_PARAMETERS_RFPARAMETERS_H_

namespace blond {
    class RfParameters;
}

#include <blond/beams/Beams.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/constants.h>
#include <blond/globals.h>
#include <blond/vector_math.h>


namespace blond {
    class RfParameters {
    public:
        enum acc_sys_t { as_single, all, first };

        int section_index;
        int &n_turns;
        double &ring_circumference;
        double &charge;
        int &alpha_order;
        f_vector_t &t_rev;
        f_vector_t &momentum;
        f_vector_t &beta;
        f_vector_t &gamma;
        f_vector_t &energy;
        f_vector_t &eta_0;
        f_vector_t &eta_1;
        f_vector_t &eta_2;


        f_vector_t E_increment;
        f_vector_t phi_s;
        f_vector_t Qs;
        f_vector_t omega_s0;
        f_vector_2d_t omega_rf_d;
        f_vector_2d_t phi_rf;
        f_vector_t dphi_rf;
        f_vector_t dphi_rf_steering;
        f_vector_t t_rf;
        f_vector_2d_t omega_rf;
        int_vector_t sign_eta_0;

        int counter;
        int n_rf;
        // Why are these 2d?
        f_vector_2d_t harmonic;
        f_vector_2d_t voltage;
        f_vector_2d_t phi_offset;
        f_vector_2d_t phi_noise;
        double length_ratio;
        double section_length;


        f_vector_t calc_phi_s(RfParameters *rfp,
                              const acc_sys_t acc_sys =
                                  acc_sys_t::as_single);


        // TODO write an eta_tracking function with a vector dE
        double eta_tracking(const Beams *beam,
                            const int counter,
                            const double dE) const;

        RfParameters(GeneralParameters *GP,
                     int _n_rf,
                     f_vector_2d_t _harmonic,
                     f_vector_2d_t _voltage,
                     f_vector_2d_t _phi_offset,
                     f_vector_2d_t _phi_noise = f_vector_2d_t(),
                     f_vector_2d_t _omega_rf = f_vector_2d_t(),
                     int _section_index = 1,
                     acc_sys_t acc_sys = as_single)
            : section_index(_section_index - 1),
              n_turns(GP->n_turns),
              ring_circumference(GP->ring_circumference),
              charge(GP->charge),
              alpha_order(GP->alpha_order),
              t_rev(GP->t_rev),
              momentum(GP->momentum[_section_index - 1]),
              beta(GP->beta[_section_index - 1]),
              gamma(GP->gamma[_section_index - 1]),
              energy(GP->energy[_section_index - 1]),
              eta_0(GP->eta_0[_section_index - 1]),
              eta_1(GP->eta_1[_section_index - 1]),
              eta_2(GP->eta_2[_section_index - 1])
        {
            counter = 0;
            n_rf = _n_rf;
            harmonic = _harmonic;
            voltage = _voltage;
            phi_noise = _phi_noise;
            phi_offset = _phi_offset;
            section_length = GP->ring_length[section_index];
            length_ratio = section_length / ring_circumference;

            for (const auto &e : eta_0)
                sign_eta_0.push_back(mymath::sign(e));

            // TODO: check with multi Rf
            E_increment.resize(n_turns);
            for (int j = 0; j < n_turns; ++j)
                E_increment[j] = energy[j + 1] - energy[j];

            phi_s = calc_phi_s(this, acc_sys);

            Qs.resize(n_turns + 1);
            for (int i = 0; i < n_turns + 1; ++i)
                Qs[i] = std::sqrt(harmonic[section_index][i] * charge
                                  * voltage[section_index][i]
                                  * std::abs(eta_0[i] * std::cos(phi_s[i])) /
                                  (2 * constant::pi * beta[i] *
                                   beta[i] * energy[i]));

            // omega_s0.resize(n_turns + 1);
            // for (int i = 0; i < (n_turns + 1); ++i)
            //     omega_s0[i] = Qs[i] * GP->omega_rev[i];
            omega_s0 = Qs * GP->omega_rev;

            // omega_rf_d.resize(n_rf, f_vector_t(n_turns + 1));
            omega_rf_d.resize(n_rf);

            for (int i = 0; i < n_rf; ++i)
                omega_rf_d[i] = (2. * constant::pi * constant::c / ring_circumference)
                                * beta * harmonic[i];
            // for (int j = 0; j < n_turns + 1; ++j)
            //     omega_rf_d[i][j] = 2. * constant::pi * beta[j] * constant::c *
            //                        harmonic[i][j] / ring_circumference;

            if (_omega_rf.empty())
                omega_rf = omega_rf_d;
            else
                omega_rf = _omega_rf;

            phi_rf = phi_offset;
            dphi_rf.resize(n_rf, 0);
            dphi_rf_steering.resize(n_rf, 0);

            t_rf = (2. * constant::pi) / omega_rf[section_index];
            // t_rf.resize(n_turns + 1);
            // for (int i = 0; i < n_turns + 1; ++i)
            //     t_rf[i] = 2 * constant::pi / omega_rf[section_index][i];
        }
        ~RfParameters() {};
    };
} // blond

#endif /* INPUT_PARAMETERS_RFPARAMETERS_H_ */
