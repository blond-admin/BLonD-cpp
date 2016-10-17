/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 14:38:41
*/

#ifndef IMPEDANCES_INDUCEDVOLTAGE_H_
#define IMPEDANCES_INDUCEDVOLTAGE_H_

#include <blond/configuration.h>
#include <blond/impedances/Intensity.h>
#include <vector>




void linear_interp_kick(const double *__restrict beam_dt,
                        double *__restrict beam_dE,
                        const double *__restrict voltage_array,
                        const double *__restrict bin_centers,
                        const int n_slices,
                        const int n_macroparticles);


class API InducedVoltage {
public:
    std::vector<double> fInducedVoltage;

    InducedVoltage() {};

    virtual void track() = 0;
    virtual void reprocess() = 0;
    virtual std::vector<double> induced_voltage_generation(uint length = 0) = 0;
    virtual ~InducedVoltage() {};
};

class API InducedVoltageTime : public InducedVoltage {
public:
    enum time_or_freq { time_domain, freq_domain };

    std::vector<Intensity *> fWakeSourceList;
    std::vector<double> fTimeArray;
    std::vector<double> fTotalWake;
    uint fCut;
    uint fShape;
    time_or_freq fTimeOrFreq;

    void track();
    void sum_wakes(std::vector<double> &v);
    void reprocess();
    std::vector<double> induced_voltage_generation(uint length = 0);
    InducedVoltageTime(std::vector<Intensity *> &WakeSourceList,
                       time_or_freq TimeOrFreq = freq_domain);

    ~InducedVoltageTime();
};

class API InducedVoltageFreq : public InducedVoltage {
public:
    enum freq_res_option_t {
        round_option,
        ceil_option,
        floor_option
    };
    // Impedance sources inputed as a list (eg: list of BBResonators objects)*
    std::vector<Intensity *> fImpedanceSourceList;

    // *Input frequency resolution in [Hz], the beam profile sampling for the
    // spectrum
    // will be adapted according to the freq_res_option.*
    double fFreqResolutionInput;

    // Number of turns to be considered as memory for induced voltage
    // calculation.*
    uint fNTurnsMem;
    bool fRecalculationImpedance;
    bool fSaveIndividualVoltages;
    // *Real frequency resolution in [Hz], according to the obtained
    // n_fft_sampling.*
    double fFreqResolution;
    // *Frequency array of the impedance in [Hz]*
    f_vector_t fFreqArray;
    uint fNFFTSampling;
    freq_res_option_t fFreqResOption;
    // *Total impedance array of all sources in* [:math:`\Omega`]
    complex_vector_t fTotalImpedance;
    complex_vector_t fMatrixSaveIndividualImpedances;
    f_vector_t fMatrixSaveIndividualVoltages;

    uint fLenArrayMem;
    uint fLenArrayMemExt;
    uint fNPointsFFT;

    f_vector_t fFreqArrayMem;
    complex_vector_t fTotalImpedanceMem;
    f_vector_t fTimeArrayMem;

    // *Induced voltage from the sum of the wake sources in [V]*
    // f_vector_t fInducedVoltage;

    void track();
    void sum_impedances(f_vector_t &);

    // Reprocess the impedance contributions with respect to the new_slicing.
    void reprocess();
    std::vector<double> induced_voltage_generation(uint length = 0);
    InducedVoltageFreq(
        std::vector<Intensity *> &impedanceSourceList,
        double freqResolutionInput = 0.0,
        freq_res_option_t freq_res_option = freq_res_option_t::round_option,
        uint NTurnsMem = 0, bool recalculationImpedance = false,
        bool saveIndividualVoltages = false);
    ~InducedVoltageFreq();
};

class API TotalInducedVoltage : public InducedVoltage {
public:
    std::vector<InducedVoltage *> fInducedVoltageList;
    std::vector<double> fTimeArray;
    std::vector<double> fRevTimeArray;
    uint fCounterTurn = 0;
    uint fNTurnsMemory;
    bool fInductiveImpedanceOn = false;

    void track();
    void track_memory();
    void track_ghosts_particles();
    std::vector<double> induced_voltage_sum(uint length = 0);
    void reprocess();

    std::vector<double> induced_voltage_generation(uint length = 0)
    {
        return std::vector<double>();
    };

    TotalInducedVoltage(std::vector<InducedVoltage *> &InducedVoltageList,
                        uint NTurnsMemory = 0,
                        std::vector<double> RevTimeArray = std::vector<double>());

    ~TotalInducedVoltage();
};

#endif /* IMPEDANCES_INDUCEDVOLTAGE_H_ */