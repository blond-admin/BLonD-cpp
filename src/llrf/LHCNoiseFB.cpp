/*
* @Author: Konstantinos Iliakis
* @Date:   2016-06-23 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-06-23 14:38:41
*
* Methods to generate RF phase noise from noise spectrum and feedback noise
* amplitude as a function of bunch length**
*
* :Authors: **Helga Timko**
*/

#include <algorithm>
#include <blond/constants.h>
#include <blond/globals.h>
#include <blond/llrf/LHCNoiseFB.h>
#include <blond/math_functions.h>

const ftype LHCNoiseFB::cfwhm = std::sqrt(2.0 / std::log(2.0));

LHCNoiseFB::LHCNoiseFB(ftype bl_target, ftype gain, ftype factor,
                       ftype update_frequency, bool variable_gain,
                       f_vector_t bunch_pattern) {
    fX = 0.0;
    fBlTarg = bl_target;
    fBlMeas = bl_target;
    fA = factor;
    fNUpdate = update_frequency;
    fVariableGain = variable_gain;
    fBunchPattern = bunch_pattern;
    auto GP = Context::GP;
    auto RfP = Context::RfP;

    if (fVariableGain) {
        fG.resize(GP->n_turns + 1);
        for (uint i = 0; i < fG.size(); ++i)
            fG[i] = gain * pow(RfP->omega_s0[0] / RfP->omega_s0[i], 2);
    } else {
        fG.resize(GP->n_turns + 1, gain);
    }

    if (fBunchPattern.empty()) {
        fFwhm = [&]() { return fwhm_single_bunch(); };
    } else {
        fFwhm = [&]() { return fwhm_multi_bunch(); };
        fBlMeasBBB.resize(fBunchPattern.size(), 0);
    }
}

LHCNoiseFB::~LHCNoiseFB() {}

// TODO test this function
void LHCNoiseFB::track() {
    auto RfP = Context::RfP;
    // Calculate PhaseNoise Feedback scaling factor as a function of measured
    // FWHM bunch length.*

    // Track only in certain turns
    if (RfP->counter % fNUpdate == 0) {
        // Update bunch length, every x turns determined in main file
        fFwhm();

        // Update noise amplitude-scaling factor
        fX = fA * fX + fG[RfP->counter] * (fBlTarg - fBlMeas);

        // Limit to range [0,1]
        fX = std::max(fX, 0.0);
        fX = std::min(fX, 1.0);

        // std::cout.precision(6);
        // std::cout << std::scientific << std::showpos;
        // std::cout << "fA " << fA << "\n";
        // std::cout << "fX " << fX << "\n";
        // std::cout << "fG[counter] " << fG[RfP->counter] << "\n";
        // std::cout << "fBlTarg " << fBlTarg << "\n";
        // std::cout << "fBlMeas " << fBlMeas << "\n";
    }
}

ftype LHCNoiseFB::fwhm_interpolation(uint_vector_t index, ftype half_height) {
    auto Slice = Context::Slice;
    const auto time_resolution = Slice->bin_centers[1] - Slice->bin_centers[0];

    const auto first = index[0];
    const auto prev = first > 0 ? first - 1 : Slice->n_slices - 1;
    const auto left =
        Slice->bin_centers[first] -
        (Slice->n_macroparticles[first] - half_height) /
            (Slice->n_macroparticles[first] - Slice->n_macroparticles[prev]) *
            time_resolution;

    const auto last = index.back();
    auto right = 0.0;
    if (last < Slice->n_slices - 1) {
        right = Slice->bin_centers[last] +
                (Slice->n_macroparticles[last] - half_height) /
                    (Slice->n_macroparticles[last] -
                     Slice->n_macroparticles[last + 1]) *
                    time_resolution;
    }

    // std::cout << "time_resolution " << time_resolution << '\n';
    // std::cout << "half_height " << half_height << '\n';
    // std::cout << "first " << first << '\n';
    // std::cout << "bin_centers[first] " << Slice->bin_centers[first] << '\n';
    // std::cout << "bin_centers[last] " << Slice->bin_centers[last] << '\n';
    // std::cout << "n_macroparticles[first] " << Slice->n_macroparticles[first]
    // << '\n';
    // std::cout << "n_macroparticles[prev] " << Slice->n_macroparticles[prev]
    // <<
    // '\n';
    // std::cout << "n_macroparticles[last] " << Slice->n_macroparticles[last]
    // <<
    // '\n';
    // std::cout << "n_macroparticles[last+1] " << Slice->n_macroparticles[last
    // +
    // 1] << '\n';
    // std::cout << "left " << left << '\n';
    // std::cout << "last " << last << '\n';
    // std::cout << "right " << right << '\n';
    // std::cout << "cfwhm " << cfwhm << '\n';

    return cfwhm * (right - left);
}

// TODO test this function
void LHCNoiseFB::fwhm_single_bunch() {
    auto Slice = Context::Slice;

    // Single-bunch FWHM bunch length calculation with interpolation.
    auto i = mymath::max(Slice->n_macroparticles.data(), Slice->n_slices);
    ftype half_height = Slice->n_macroparticles[i] / 2;

    uint_vector_t index;

    for (uint i = 0; i < Slice->n_slices; ++i)
        if (Slice->n_macroparticles[i] > half_height)
            index.push_back(i);

    if (index.size() == 0) {
        // std::cerr << "[LHCNoiseFB] ERROR! index vector should have at least
        // one element\n";
        return;
    }

    fBlMeas = fwhm_interpolation(index, half_height);
}

// TODO test this function
void LHCNoiseFB::fwhm_multi_bunch() {

    // Multi-bunch FWHM bunch length calculation with interpolation.*

    auto Slice = Context::Slice;
    auto RfP = Context::RfP;

    // Find correct RF buckets
    auto phi_RF = RfP->phi_RF[0][RfP->counter];
    auto omega_RF = RfP->omega_RF[0][RfP->counter];
    // std::cout << phi_RF << "\n";
    // std::cout << omega_RF << "\n";
    f_vector_t bucket_min(fBunchPattern.size());
    for (uint i = 0; i < bucket_min.size(); ++i)
        bucket_min[i] =
            (phi_RF + 2 * constant::pi * fBunchPattern[i]) / omega_RF;

    f_vector_t bucket_max(bucket_min.size());
    for (uint i = 0; i < bucket_min.size(); ++i)
        bucket_max[i] = bucket_min[i] + 2 * constant::pi / omega_RF;

    // Bunch-by-bunch FWHM bunch length
    for (uint i = 0; i < fBunchPattern.size(); ++i) {
        uint_vector_t bind;
        for (uint j = 0; j < Slice->n_slices; ++j) {
            auto val = (Slice->bin_centers[j] - bucket_min[i]) *
                           (Slice->bin_centers[j] - bucket_max[i]) <
                       0;
            if (val)
                bind.push_back(j);
        }

        ftype hheight = 0;
        for (const auto& j : bind) {
            if (Slice->n_macroparticles[j] > hheight)
                hheight = Slice->n_macroparticles[j];
        }
        // std::cout << "height: " << hheight << "\n";
        uint_vector_t index;

        uint k = 0;
        for (const auto& j : bind) {
            if (Slice->n_macroparticles[j] > hheight)
                index.push_back(bind[k]);
            k++;
        }
        // std::cout << "index size = " << index.size() << "\n";
        if (index.empty()) {
            // std::cerr << "[LHCNoiseFB] ERROR! index vector should have at
            // least one element\n";
            continue;
        }
        fBlMeasBBB[i] = fwhm_interpolation(index, hheight);
    }

    fBlMeas = mymath::mean(fBlMeasBBB.data(), fBlMeasBBB.size());
}
