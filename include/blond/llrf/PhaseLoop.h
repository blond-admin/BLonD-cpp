/*
 * PhaseLoop.h
 *
 *  Created on: Apr 7, 2016
 *      Author: kiliakis
 */

#ifndef PHASELOOP_H_
#define PHASELOOP_H_
#include <blond/utilities.h>
#include <blond/configuration.h>
namespace blond {

	class API PhaseLoop {
	public:
		virtual void track() {};
		void default_track();
		PhaseLoop(ftype *PL_gain, ftype window_coefficient, int _delay,
			ftype *_phaseNoise, ftype *_LHCNoiseFB);
		void beam_phase();
		void phase_difference();
		void radial_steering_from_freq();
		PhaseLoop() {};
		int delay = 0;
		ftype alpha = 0;
		ftype *gain = NULL;
		ftype drho = 0;
		ftype domega_RF = 0;
		ftype phi_beam = 0;
		ftype dphi = 0;
		ftype reference = 0;
		ftype *RFnoise = NULL;
		ftype *noiseFB = NULL;
		virtual ~PhaseLoop() {};
	};

	class API LHC : public PhaseLoop {
	private:
		//ftype gain;
		//ftype domega_RF;
	public:
		ftype gain2;
		ftype lhc_y;
		ftype *lhc_a;
		ftype *lhc_t;

		~LHC();
		void track();
		LHC(ftype *PL_gain, ftype SL_gain = 0, ftype window_coefficient = 0,
			ftype *phaseNoise = NULL, ftype *LHCNoiseFB = NULL, int delay = 0);
	};

	class API PSB : public PhaseLoop {
	private:
		ftype *gain2;
		//ftype *gain;
		int *dt;
		ftype average_dE;
		int PL_counter;
		std::vector<int> on_time;
		ftype *coefficients;
		ftype dphi_av;
		ftype dphi_av_prev;
		ftype drho_prev;
		ftype t_accum;
		ftype domega_PL;
		ftype domega_RL;
		//ftype domega_RF;
	public:
		~PSB();
		void track();
		PSB(ftype *PL_gain, ftype *RL_gain = NULL, ftype PL_period = 0,
			ftype RL_period = 0, ftype *coefficients = NULL,
			ftype window_coefficient = 0, ftype *phaseNoise = NULL,
			ftype *LHCNoiseFB = NULL, int delay = 0);
		void radial_difference();
		void precalculate_time();
	};
}
#endif /* PHASELOOP_H_ */
