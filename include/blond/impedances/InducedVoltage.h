/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 14:38:41
*/

#ifndef IMPEDANCES_INDUCEDVOLTAGE_H_
#define IMPEDANCES_INDUCEDVOLTAGE_H_

#include <vector>
#include <blond/utilities.h>
#include <blond/configuration.h>
#include <blond/impedances/Intensity.h>

//#include <complex>
//typedef std::complex<float> complex_t;
namespace blond {

	enum time_or_freq {
		time_domain, freq_domain
	};
	//unsigned int next_regular(unsigned int target);



	class API InducedVoltage {
	public:
		std::vector<ftype> fInducedVoltage;

		InducedVoltage() {};
		inline void linear_interp_kick(const ftype *__restrict beam_dt,
			ftype *__restrict beam_dE,
			const ftype *__restrict voltage_array,
			const ftype *__restrict bin_centers,
			const int n_slices,
			const int n_macroparticles,
			const ftype acc_kick = 0.0);
		virtual void track() = 0;
		virtual void reprocess() = 0;
		virtual std::vector<ftype> induced_voltage_generation(unsigned int length = 0) = 0;
		virtual ~InducedVoltage() {};
	};



	class API InducedVoltageTime : public InducedVoltage {
	public:


		std::vector<Intensity *> fWakeSourceList;
		std::vector<ftype> fTimeArray;
		std::vector<ftype> fTotalWake;
		unsigned int fCut;
		unsigned int fShape;
		time_or_freq fTimeOrFreq;


		void track();
		void sum_wakes(std::vector<ftype> &v);
		void reprocess();
		std::vector<ftype> induced_voltage_generation(unsigned int length = 0);
		InducedVoltageTime(std::vector<Intensity *> &WakeSourceList,
			time_or_freq TimeOrFreq = freq_domain);
		~InducedVoltageTime() {};
	};


	class InducedVoltageFreq : public InducedVoltage {
	public:
		void track();
		void sum_impedances();
		void reprocess();
		std::vector<ftype> induced_voltage_generation(unsigned int length = 0);
		InducedVoltageFreq();
		~InducedVoltageFreq() {};
	};


	class API TotalInducedVoltage : public InducedVoltage {
	public:
		std::vector<InducedVoltage *> fInducedVoltageList;
		std::vector<ftype> fTimeArray;
		std::vector<ftype> fRevTimeArray;
		unsigned int fCounterTurn = 0;
		unsigned int fNTurnsMemory;
		bool fInductiveImpedanceOn = false;

		void track();
		void track_memory();
		void track_ghosts_particles();
		std::vector<ftype> induced_voltage_sum(unsigned int length = 0);
		void reprocess();
		std::vector<ftype> induced_voltage_generation(unsigned int length = 0);

		TotalInducedVoltage(std::vector<InducedVoltage *> &InducedVoltageList,
			unsigned int NTurnsMemory = 0,
			std::vector<ftype> RevTimeArray = std::vector<ftype>());
		~TotalInducedVoltage() {};

	};
}
#endif /* IMPEDANCES_INDUCEDVOLTAGE_H_ */