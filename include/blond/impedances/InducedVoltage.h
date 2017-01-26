/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 14:38:41
*/

#ifndef IMPEDANCES_INDUCEDVOLTAGE_H_
#define IMPEDANCES_INDUCEDVOLTAGE_H_

#include <blond/configuration.h>
#include <blond/beams/Beams.h>
#include <blond/impedances/Intensity.h>
#include <vector>

namespace blond {
    void linear_interp_kick(const double *__restrict beam_dt,
                            double *__restrict beam_dE,
                            const double *__restrict voltage_array,
                            const double *__restrict bin_centers,
                            const int n_slices,
                            const int n_macroparticles);


    class InducedVoltage {
    public:
        f_vector_t fInducedVoltage;
        Slices *fSlices;

        InducedVoltage() {};

        virtual void track(Beams *beam) = 0;
        virtual void reprocess(Slices *slices) = 0;
        virtual f_vector_t induced_voltage_generation(Beams *beam,
                int length = 0) = 0;
        virtual ~InducedVoltage() {};
    };

    class InducedVoltageTime : public InducedVoltage {
    public:
        enum time_or_freq { time_domain, freq_domain };

        std::vector<Intensity *> fWakeSourceList;
        f_vector_t fTimeArray;
        f_vector_t fTotalWake;
        int fCut;
        int fShape;
        time_or_freq fTimeOrFreq;

        void track(Beams *beam);
        void sum_wakes(f_vector_t &v);
        void reprocess(Slices *newSlices);
        f_vector_t induced_voltage_generation(Beams *beam, int length = 0);
        InducedVoltageTime(Slices *slices,
                           const std::vector<Intensity *> &WakeSourceList,
                           time_or_freq TimeOrFreq = freq_domain);

        ~InducedVoltageTime();
    };

    class InducedVoltageFreq : public InducedVoltage {
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
        int fNTurnsMem;
        bool fRecalculationImpedance;
        bool fSaveIndividualVoltages;
        // *Real frequency resolution in [Hz], according to the obtained
        // n_fft_sampling.*
        double fFreqResolution;
        // *Frequency array of the impedance in [Hz]*
        f_vector_t fFreqArray;
        int fNFFTSampling;
        freq_res_option_t fFreqResOption;
        // *Total impedance array of all sources in* [:math:`\Omega`]
        complex_vector_t fTotalImpedance;
        complex_vector_t fMatrixSaveIndividualImpedances;
        f_vector_t fMatrixSaveIndividualVoltages;

        int fLenArrayMem;
        int fLenArrayMemExt;
        int fNPointsFFT;

        f_vector_t fFreqArrayMem;
        complex_vector_t fTotalImpedanceMem;
        f_vector_t fTimeArrayMem;

        // *Induced voltage from the sum of the wake sources in [V]*
        // f_vector_t fInducedVoltage;

        void track(Beams *beam);
        void sum_impedances(f_vector_t &);

        // Reprocess the impedance contributions with respect to the new_slicing.
        void reprocess(Slices *newSlices);
        f_vector_t induced_voltage_generation(Beams *beam, int length = 0);
        InducedVoltageFreq(Slices *slices,
                           const std::vector<Intensity *> &impedanceSourceList,
                           double freqResolutionInput = 0.0,
                           freq_res_option_t freq_res_option = freq_res_option_t::round_option,
                           int NTurnsMem = 0, bool recalculationImpedance = false,
                           bool saveIndividualVoltages = false);
        ~InducedVoltageFreq();
    };

    class TotalInducedVoltage : public InducedVoltage {
    public:
        std::vector<InducedVoltage *> fInducedVoltageList;
        f_vector_t fTimeArray;
        f_vector_t fRevTimeArray;
        int fCounterTurn = 0;
        int fNTurnsMemory;
        bool fInductiveImpedanceOn = false;

        Beams *fBeam;

        void track(Beams *beam);
        void track_memory();
        void track_ghosts_particles(Beams *ghostBeam);
        f_vector_t induced_voltage_sum(Beams *beam, int length = 0);
        void reprocess(Slices *newSlices);

        f_vector_t induced_voltage_generation(Beams *beam, int length = 0)
        {
            return f_vector_t();
        };

        TotalInducedVoltage(Beams *beam,
                            Slices *slices,
                            const std::vector<InducedVoltage *> &InducedVoltageList,
                            int NTurnsMemory = 0,
                            f_vector_t RevTimeArray = f_vector_t());

        ~TotalInducedVoltage();
    };
} // blond
#endif /* IMPEDANCES_INDUCEDVOLTAGE_H_ */