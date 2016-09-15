/*
 * RfParameters.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef INPUT_PARAMETERS_RFPARAMETERS_H_
#define INPUT_PARAMETERS_RFPARAMETERS_H_


class RfParameters;

#include <blond/beams/Beams.h>
#include <blond/input_parameters/GeneralParameters.h>


class API RfParameters {
public:
    enum accelerating_systems_t { as_single, all, first };

    f_vector_t E_increment;
    f_vector_t phi_s;
    f_vector_t Qs;
    f_vector_t omega_s0;
    f_vector_2d_t omega_RF_d;
    f_vector_2d_t phi_RF;
    f_vector_t dphi_RF;
    f_vector_t dphi_RF_steering;
    f_vector_t t_RF;
    f_vector_2d_t omega_RF;

    // TODO assume input_value is an array
    // that is why we don't have any input_check function
    uint counter;
    uint n_rf;
    // int n_turns;
    f_vector_2d_t harmonic;
    f_vector_2d_t voltage;
    f_vector_2d_t phi_offset;
    f_vector_2d_t phi_noise;
    uint idx;
    ftype length_ratio;
    ftype section_length;


    // TODO write an eta_tracking function with a vector dE
    ftype eta_tracking(const Beams *beam,
                       const uint counter,
                       const ftype dE) const;
    ftype eta_0(const uint i) const;
    ftype eta_1(const uint i) const;
    ftype eta_2(const uint i) const;
    ftype beta(const uint i) const;
    ftype gamma(const uint i) const;
    ftype energy(const uint i) const;
    ftype momentum(const uint i) const;
    int sign_eta_0(const uint i) const;

    RfParameters(uint _n_rf, f_vector_2d_t _harmonic, f_vector_2d_t _voltage,
                 f_vector_2d_t _phi_offset,
                 f_vector_2d_t _phi_noise = f_vector_2d_t(),
                 f_vector_2d_t _omega_rf = f_vector_2d_t(),
                 uint _section_index = 1,
                 accelerating_systems_t accelerating_systems = as_single);
    ~RfParameters();

private:
};

void calc_phi_s(ftype *out, RfParameters *rfp,
                const RfParameters::accelerating_systems_t acc_sys =
                    RfParameters::accelerating_systems_t::as_single);

#endif /* INPUT_PARAMETERS_RFPARAMETERS_H_ */
