
#include <blond/constants.h>
#include <blond/globals.h>
#include <blond/impedances/InducedVoltage.h>
#include <blond/math_functions.h>
#include <blond/fft.h>
#include <blond/utilities.h>
#include <blond/vector_math.h>

void linear_interp_kick(
    const double *__restrict beam_dt,
    double *__restrict beam_dE,
    const double *__restrict voltage_array,
    const double *__restrict bin_centers,
    const int n_slices,
    const int n_macroparticles)
{

    const double binFirst = bin_centers[0];
    const double binLast = bin_centers[n_slices - 1];
    const double inv_bin_width = (n_slices - 1) / (binLast - binFirst);

    #pragma omp parallel for
    for (int i = 0; i < n_macroparticles; ++i) {
        const double a = beam_dt[i];
        const int ffbin = static_cast<int>((a - binFirst) * inv_bin_width);
        const double voltageKick =
            ((a < binFirst) || (a > binLast))
            ? 0.0
            : voltage_array[ffbin] +
            (a - bin_centers[ffbin]) *
            (voltage_array[ffbin + 1] - voltage_array[ffbin]) *
            inv_bin_width;
        beam_dE[i] += voltageKick;
    }
}

InducedVoltageTime::InducedVoltageTime(Slices *slices,
                                       const std::vector<Intensity *> &WakeList,
                                       time_or_freq TimeOrFreq)
{
    // Induced voltage derived from the sum of
    // several wake fields (time domain).*
    fSlices = slices;
    // *Wake sources inputed as a list (eg: list of BBResonators objects)*
    fWakeSourceList = WakeList;

    // *Time array of the wake in [s]*
    fTimeArray = f_vector_t();

    // *Total wake array of all sources in* [:math:`\Omega / s`]
    fTotalWake = f_vector_t();

    // *Induced voltage from the sum of the wake sources in [V]*
    fInducedVoltage = f_vector_t();

    // Pre-processing the wakes
    fTimeArray = fSlices->bin_centers - fSlices->bin_centers[0];
    sum_wakes(fTimeArray);

    fCut = fTimeArray.size() + fSlices->n_slices - 1;
    fShape = mymath::next_regular(fCut);

    fTimeOrFreq = TimeOrFreq;
}

InducedVoltageTime::~InducedVoltageTime() { fft::destroy_plans(); }

inline void InducedVoltageTime::track(Beams *beam)
{
    // Tracking Method
    f_vector_t v = this->induced_voltage_generation(beam) * beam->charge;

    linear_interp_kick(beam->dt.data(), beam->dE.data(), v.data(),
                       fSlices->bin_centers.data(), fSlices->n_slices,
                       beam->n_macroparticles);
}

void InducedVoltageTime::sum_wakes(f_vector_t &TimeArray)
{
    // *Summing all the wake contributions in one total wake.*
    fTotalWake.resize(TimeArray.size());
    std::fill(fTotalWake.begin(), fTotalWake.end(), 0);
    for (const auto &i : fWakeSourceList) {
        i->wake_calc(TimeArray);
        fTotalWake += i->fWake;
    }
}

void InducedVoltageTime::reprocess(Slices *newSlices)
{
    // *Reprocess the wake contributions with respect to the new_slicing.*
    // WARNING As Slice is a global variable,
    // users will have to change this variable and call reprocess()

    fSlices = newSlices;
    // fTimeArray.resize(fSlices->n_slices);
    fTimeArray = fSlices->bin_centers - fSlices->bin_centers[0];
    sum_wakes(fTimeArray);

    fCut = fTimeArray.size() + fSlices->n_slices - 1;
    fShape = mymath::next_regular(fCut);
}

f_vector_t InducedVoltageTime::induced_voltage_generation(Beams *beam,
        uint length)
{

    // Method to calculate the induced voltage from wakes with convolution.*
    f_vector_t inducedVoltage;

    const double factor = -beam->charge * constant::e * beam->intensity
                          / beam->n_macroparticles;

    if (fTimeOrFreq == freq_domain) {
        auto in1 = fSlices->n_macroparticles;
        auto in2 = fTotalWake;

        fft::convolution_with_ffts(in1, in2, inducedVoltage);
        inducedVoltage *= factor;

    } else if (fTimeOrFreq == time_domain) {

        inducedVoltage.resize(fTotalWake.size() + fSlices->n_slices - 1);

        mymath::convolution(fTotalWake.data(), fTotalWake.size(),
                            fSlices->n_macroparticles.data(),
                            fSlices->n_macroparticles.size(),
                            inducedVoltage.data());
        inducedVoltage *= factor;

    } else {
        std::cerr << "Error: Only freq_domain or time_domain are allowed\n";
        exit(-1);
    }

    fInducedVoltage = inducedVoltage;
    fInducedVoltage.resize((uint)fSlices->n_slices);

    if (length > 0)
        inducedVoltage.resize(length, 0);

    return inducedVoltage;
}

InducedVoltageFreq::InducedVoltageFreq(Slices *slices,
                                       const std::vector<Intensity *> &impedList,
                                       double freqResolutionInput,
                                       freq_res_option_t freq_res_option,
                                       uint NTurnsMem,
                                       bool recalculationImpedance,
                                       bool saveIndividualVoltages)
{
    fNTurnsMem = NTurnsMem;
    fSlices = slices;
    fImpedanceSourceList = impedList;
    fFreqResolutionInput = freqResolutionInput;

    // *Length of one slice.*
    auto timeResolution = (fSlices->bin_centers[1] - fSlices->bin_centers[0]);
    fRecalculationImpedance = recalculationImpedance;
    fFreqResOption = freq_res_option;

    if (fNTurnsMem == 0) {

        if (fFreqResolutionInput == 0) {
            fNFFTSampling = fSlices->n_slices;
        } else {
            int a;
            double b = 1 / (fFreqResolutionInput * timeResolution);
            switch (fFreqResOption) {
                case freq_res_option_t::round_option:
                    a = std::round(b);
                    break;
                case freq_res_option_t::ceil_option:
                    a = std::ceil(b);
                    break;
                case freq_res_option_t::floor_option:
                    a = std::floor(b);
                    break;
                default:
                    std::cerr << "The input freq_res_option is not recognized\n";
                    exit(-1);
                    break;
            }
            fNFFTSampling = mymath::next_regular(a);

            if ((int) fNFFTSampling < fSlices->n_slices) {
                std::cerr << "The input frequency resolution step is too big, "
                          "and the whole\n"
                          << "bunch is not sliced... The number of sampling "
                          "points for the\n"
                          << "FFT is corrected in order to sample the whole "
                          "bunch (and\n"
                          << "you might consider changing the input in order "
                          "to have\n"
                          << "a finer resolution\n";
                fNFFTSampling = mymath::next_regular(fSlices->n_slices);
            }
        }

        fFreqResolution = 1 / (fNFFTSampling * timeResolution);

        fFreqArray = fft::rfftfreq(fNFFTSampling, timeResolution);
        sum_impedances(fFreqArray);

        fSaveIndividualVoltages = saveIndividualVoltages;
        if (fSaveIndividualVoltages) {
            // Do I really need to store the length??
            uint length = fImpedanceSourceList.size();
            fMatrixSaveIndividualImpedances =
                complex_vector_t(length * fFreqArray.size(), 0);
            fMatrixSaveIndividualVoltages =
                f_vector_t(length * fSlices->n_slices, 0);
            for (uint i = 0; i < length; ++i) {
                const uint row_width =
                    fImpedanceSourceList[i]->fImpedance.size();
                for (uint j = 0; j < row_width; ++j) {
                    fMatrixSaveIndividualImpedances[i * row_width + j] =
                        fImpedanceSourceList[i]->fImpedance[j];
                }
            }
        }

    } else {
        fNTurnsMem = NTurnsMem;
        fLenArrayMem = (fNTurnsMem + 1) * fSlices->n_slices;
        fLenArrayMemExt = (fNTurnsMem + 2) * fSlices->n_slices;
        fNPointsFFT = mymath::next_regular(fLenArrayMemExt);
        fFreqArrayMem = fft::rfftfreq(fNPointsFFT, timeResolution);
        fTotalImpedanceMem =
            complex_vector_t(fFreqArrayMem.size(), complex_t(0, 0));

        fTimeArrayMem.reserve((fNTurnsMem + 1) * fSlices->n_slices);
        const double factor = fSlices->edges.back() - fSlices->edges.front();

        for (uint i = 0; i < fNTurnsMem + 1; ++i) {
            for (int j = 0; j < fSlices->n_slices; ++j) {
                fTimeArrayMem.push_back(fSlices->bin_centers[j] + factor * i);
            }
        }

        for (const auto &impObj : fImpedanceSourceList) {
            impObj->imped_calc(fFreqArrayMem);
            fTotalImpedanceMem += impObj->fImpedance;
        }
    }
}

InducedVoltageFreq::~InducedVoltageFreq() { fft::destroy_plans(); }

void InducedVoltageFreq::track(Beams *beam)
{
    // Tracking Method

    induced_voltage_generation(beam);
    auto v = fInducedVoltage * beam->charge;

    linear_interp_kick(beam->dt.data(), beam->dE.data(), v.data(),
                       fSlices->bin_centers.data(), fSlices->n_slices,
                       beam->n_macroparticles);
}

void InducedVoltageFreq::sum_impedances(f_vector_t &freq_array)
{

    fTotalImpedance.resize(freq_array.size());
    std::fill(fTotalImpedance.begin(), fTotalImpedance.end(), complex_t(0, 0));
    for (const auto &i : fImpedanceSourceList) {
        i->imped_calc(freq_array);
        fTotalImpedance += i->fImpedance;
    }
}

void InducedVoltageFreq::reprocess(Slices *newSlices)
{
    fSlices = newSlices;

    auto timeResolution = (fSlices->bin_centers[1] - fSlices->bin_centers[0]);
    if (fFreqResolutionInput == 0) {
        fNFFTSampling = fSlices->n_slices;
    } else {
        int a;
        double b = 1 / (fFreqResolutionInput * timeResolution);
        switch (fFreqResOption) {
            case freq_res_option_t::round_option:
                a = std::round(b);
                break;
            case freq_res_option_t::ceil_option:
                a = std::ceil(b);
                break;
            case freq_res_option_t::floor_option:
                a = std::floor(b);
                break;
            default:
                std::cerr << "The input freq_res_option is not recognized\n";
                exit(-1);
                break;
        }
        fNFFTSampling = mymath::next_regular(a);

        if ((int) fNFFTSampling < fSlices->n_slices) {
            std::cerr
                    << "The input frequency resolution step is too big, and the "
                    "whole\n"
                    << "bunch is not sliced... The number of sampling points for "
                    "the\n"
                    << "FFT is corrected in order to sample the whole bunch (and\n"
                    << "you might consider changing the input in order to have\n"
                    << "a finer resolution\n";
            fNFFTSampling = mymath::next_regular(fSlices->n_slices);
        }
    }

    fFreqResolution = 1 / (fNFFTSampling * timeResolution);

    fSlices->beam_spectrum_generation(fNFFTSampling, true);
    fFreqArray = fSlices->fBeamSpectrumFreq;

    fTotalImpedance.clear();
    sum_impedances(fFreqArray);
}

f_vector_t InducedVoltageFreq::induced_voltage_generation(Beams *beam,
        uint length)
{
    //    Method to calculate the induced voltage from the inverse FFT of the
    //    impedance times the spectrum (fourier convolution).

    if (fRecalculationImpedance)
        sum_impedances(fFreqArray);

    fSlices->beam_spectrum_generation(fNFFTSampling);
    const auto n = fImpedanceSourceList.size();
    const auto factor = -beam->charge * constant::e * beam->ratio *
                        fSlices->fBeamSpectrumFreq[1] * 2 *
                        (fSlices->fBeamSpectrum.size() - 1);

    if (fSaveIndividualVoltages) {

        for (uint i = 0; i < n; ++i) {
            f_vector_t res;
            complex_vector_t in(fSlices->fBeamSpectrum.size());

            for (uint j = 0; j < in.size(); ++j) {
                in[j] = fMatrixSaveIndividualImpedances[j * n + i] *
                        fSlices->fBeamSpectrum[j];
            }

            fft::irfft(in, res, 0, Context::n_threads);

            assert((int)res.size() >= fSlices->n_slices);

            res.resize(fSlices->n_slices);
            res *= factor;

            for (int j = 0; j < fSlices->n_slices; ++j) {
                fMatrixSaveIndividualVoltages[j * n + i] = res[j];
            }
        }

        fInducedVoltage.clear();
        fInducedVoltage.resize(fSlices->fBeamSpectrum.size());
        for (uint i = 0; i < fSlices->fBeamSpectrum.size(); ++i) {
            double sum = 0.0;
            for (uint j = 0; j < n; ++j) {
                sum += fMatrixSaveIndividualVoltages[i * n + j];
            }
            fInducedVoltage[i] = sum;
        }

        return f_vector_t();

    } else {
        f_vector_t res;
        complex_vector_t in(fSlices->fBeamSpectrum.size());
        for (uint j = 0; j < in.size(); ++j) {
            in[j] = fTotalImpedance[j] * fSlices->fBeamSpectrum[j];
        }

        fft::irfft(in, res, 0, Context::n_threads);
        assert((int)res.size() >= fSlices->n_slices);

        res.resize(fSlices->n_slices);
        res *= factor;

        fInducedVoltage = res;

        if (length > 0) {
            if (length > res.size())
                res.resize(length, 0);
            else
                res.resize(length);
        }

        return res;
    }
}

TotalInducedVoltage::TotalInducedVoltage(Beams *beam, Slices *slices,
        const std::vector<InducedVoltage *> &InducedVoltageList,
        uint NTurnsMemory,
        f_vector_t RevTimeArray)
{
    fBeam = beam;
    fSlices = slices;
    fInducedVoltageList = InducedVoltageList;
    fNTurnsMemory = NTurnsMemory;
    fInducedVoltage = f_vector_t();
    fTimeArray = fSlices->bin_centers;
}

TotalInducedVoltage::~TotalInducedVoltage() { fft::destroy_plans(); }

void TotalInducedVoltage::track(Beams *beam)
{

    this->induced_voltage_sum(beam);
    auto v = this->fInducedVoltage * beam->charge;

    linear_interp_kick(beam->dt.data(), beam->dE.data(), v.data(),
                       fSlices->bin_centers.data(), fSlices->n_slices,
                       beam->n_macroparticles);
}

void TotalInducedVoltage::track_memory() {}

void TotalInducedVoltage::track_ghosts_particles(Beams *ghostBeam) {}

void TotalInducedVoltage::reprocess(Slices *newSlices)
{
    fSlices = newSlices;

    for (auto &v : fInducedVoltageList)
        v->reprocess(newSlices);
}

f_vector_t TotalInducedVoltage::induced_voltage_sum(Beams *beam, uint length)
{
    // Method to sum all the induced voltages in one single array.
    f_vector_t tempIndVolt;
    f_vector_t extIndVolt;

    for (auto &v : fInducedVoltageList) {
        auto a = v->induced_voltage_generation(beam, length);

        if (length > 0) {
            extIndVolt.resize(a.size(), 0);
            extIndVolt += a;
        }
        tempIndVolt.resize(v->fInducedVoltage.size(), 0);
        tempIndVolt += v->fInducedVoltage;
    }

    fInducedVoltage = tempIndVolt;
    return extIndVolt;
}
