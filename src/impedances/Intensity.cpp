/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 15:25:33
*/

#include <blond/constants.h>
#include <blond/impedances/Intensity.h>
#include <blond/math_functions.h>
#include <blond/vector_math.h>

Resonators::Resonators(f_vector_t &RS, f_vector_t &FrequencyR, f_vector_t &Q)
{
    fRS = RS;
    fFrequencyR = FrequencyR;
    fQ = Q;
    fNResonators = RS.size();
    fOmegaR = 2. * constant::pi * fFrequencyR;
}

Resonators::~Resonators() {}

void Resonators::wake_calc(const f_vector_t &NewTimeArray)
{
    /*
    * Wake calculation method as a function of time.*
    */
    fTimeArray = NewTimeArray;
    fWake.resize(fTimeArray.size());
    std::fill_n(fWake.begin(), fWake.size(), 0);

    for (uint i = 0; i < fNResonators; ++i) {
        double alpha = fOmegaR[i] / (2 * fQ[i]);
        double omega_bar = std::sqrt(fOmegaR[i] * fOmegaR[i] - alpha * alpha);

        for (uint j = 0; j < fWake.size(); ++j) {
            double temp = fTimeArray[j];
            // util::dump(&temp, 1, "temp ");
            int sign = (temp > 0) - (temp < 0);
            fWake[j] += (sign + 1) * fRS[i] * alpha * std::exp(-alpha * temp) *
                        (std::cos(omega_bar * temp) -
                         alpha / omega_bar * std::sin(omega_bar * temp));
        }
    }
}

void Resonators::imped_calc(const f_vector_t &NewFrequencyArray)
{
    /*
    * Impedance calculation method as a function of frequency.*
    */

    fFreqArray = NewFrequencyArray;
    fImpedance.resize(fFreqArray.size());
    std::fill_n(fImpedance.begin(), fImpedance.size(), complex_t(0, 0));

    for (uint i = 0; i < fNResonators; ++i) {
        for (uint j = 1; j < fImpedance.size(); ++j) {
            fImpedance[j] +=
                complex_t(fRS[i], 0) /
                complex_t(1, fQ[i] * (fFreqArray[j] / fFrequencyR[i] -
                                      fFrequencyR[i] / fFreqArray[j]));
        }
    }
}

InputTable::InputTable(const f_vector_t &input1, const f_vector_t &input2,
                       const f_vector_t input3)
{
    if (input3.empty()) {

        fTimeArray = input1;
        fWakeArray = input2;

    } else {
        fFrequencyArrayLoaded = input1;
        fReZArrayLoaded = input2;
        fImZArrayLoaded = input3;
        assert(fReZArrayLoaded.size() == fImZArrayLoaded.size());

        for (uint i = 0; i < fReZArrayLoaded.size(); ++i) {
            complex_t z(fReZArrayLoaded[i], fImZArrayLoaded[i]);
            fImpedanceLoaded.push_back(z);
        }

        if (fFrequencyArrayLoaded[0] != 0) {

            fFrequencyArrayLoaded.insert(fFrequencyArrayLoaded.begin(), 0);
            fReZArrayLoaded.insert(fReZArrayLoaded.begin(), 0);
            fImZArrayLoaded.insert(fImZArrayLoaded.begin(), 0);
        }
    }
}

InputTable::~InputTable() {}

void InputTable::wake_calc(const f_vector_t &NewTimeArray)
{
    mymath::interp(NewTimeArray, fTimeArray, fWakeArray, fWake, 0.0f, 0.0f);
}

void InputTable::imped_calc(const f_vector_t &NewFrequencyArray)
{
    // Impedance calculation method as a function of frequency.*
    f_vector_t ReZ;
    f_vector_t ImZ;

    mymath::interp(NewFrequencyArray, fFrequencyArrayLoaded,
                   fReZArrayLoaded, ReZ, 0.0f, 0.0f);

    mymath::interp(NewFrequencyArray, fFrequencyArrayLoaded,
                   fImZArrayLoaded, ImZ, 0.0f, 0.0f);

    fFreqArray = NewFrequencyArray;

    // Initializing real and imaginary part separately has been
    // omitted
    fImpedance.resize(ReZ.size());
    for (uint i = 0; i < ReZ.size(); ++i) {
        fImpedance[i] = complex_t(ReZ[i], ImZ[i]);
    }
}



TravelingWaveCavity::TravelingWaveCavity(f_vector_t &RS, f_vector_t &FrequencyR, f_vector_t &aFactor)
{
    fRS = RS;
    fFrequencyR = FrequencyR;
    faFactor = aFactor;
    fNResonators = RS.size();
    // fOmegaR = 2. * constant::pi * fFrequencyR;
}

TravelingWaveCavity::~TravelingWaveCavity() {}

void TravelingWaveCavity::wake_calc(const f_vector_t &NewTimeArray)
{
    /*
     * Wake calculation method as a function of time.*
     */

    fTimeArray = NewTimeArray;
    fWake.resize(fTimeArray.size());
    std::fill_n(fWake.begin(), fWake.size(), 0);

    for (uint i = 0; i < fNResonators; ++i) {
        double a_tilde = faFactor[i] / (2 * constant::pi);
        for (uint j = 0; j < fWake.size(); ++j) {
            if (fTimeArray[j] <= a_tilde) {
                fWake[j] += (mymath::sign(fTimeArray[j]) + 1) *
                            2. * fRS[i] / a_tilde *
                            (1. - fTimeArray[j] / a_tilde) *
                            std::cos(2. * constant::pi * fFrequencyR[i]
                                     * fTimeArray[j]);
            }
        }
    }
}

void TravelingWaveCavity::imped_calc(const f_vector_t &NewFrequencyArray)
{
    /*
     * Impedance calculation method as a function of frequency.*
     */

    // TODO, make unit-tests

    fFreqArray = NewFrequencyArray;
    fImpedance.resize(fFreqArray.size());
    std::fill_n(fImpedance.begin(), fImpedance.size(), complex_t(0, 0));

    for (uint i = 0; i < fNResonators; i++) {
        for (uint j = 0; j < fImpedance.size(); j++) {
            auto minusFactor = faFactor[i] * (fFreqArray[j] - fFrequencyR[i]);
            auto plusFactor = faFactor[i] * (fFreqArray[j] + fFrequencyR[i]);

            auto Zplus = fRS[i] * complex_t(
                             pow(std::sin(minusFactor / 2.)
                                 / (minusFactor / 2.), 2),
                             -2. * (minusFactor - std::sin(minusFactor)) /
                             pow(minusFactor, 2)
                         );
            auto Zminus = fRS[i] * complex_t(
                              pow(std::sin(plusFactor / 2.)
                                  / (plusFactor / 2.), 2),
                              -2. * (plusFactor - std::sin(plusFactor)) /
                              pow(plusFactor, 2)
                          );
            fImpedance[j] += Zplus + Zminus;
        }
    }
}
