/*
 * Beams.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#include "../input_parameters/GeneralParameters.h"

#include "../includes/utilities.h"

#ifndef BEAMS_BEAMS_H_
#define BEAMS_BEAMS_H_

//#include "../includes/globals.h"

class Beams {
public:
	ftype *dt;
	ftype *dE;
	ftype mean_dt;
	ftype mean_dE;
	ftype sigma_dt;
	ftype sigma_dE;
	ftype ratio;
	ftype epsn_rms_l;
	int n_macroparticles_lost;
	int n_macroparticles;

	GeneralParameters *gp;
	int *id;
	Beams(GeneralParameters *gp, int _n_macroparticles, int _intensity);
	int n_macroparticles_alive();
	void losses_longitudinal_cut(ftype dt_min, ftype dt_max);
	void losses_energy_cut(ftype dE_min, ftype dE_max);

private:
	int intensity;
	void statistics();
};

Beams::Beams(GeneralParameters *_gp, int _n_macroparticles, int _intensity) {
	// TODO many variables are never used
	this->gp = _gp;
	this->n_macroparticles = _n_macroparticles;
	this->intensity = _intensity;
	this->dt = new ftype[n_macroparticles];
	this->dE = new ftype[n_macroparticles];
	this->mean_dt = this->mean_dE = 0;
	this->sigma_dt = this->sigma_dE = 0;
	this->ratio = intensity / n_macroparticles;
	this->epsn_rms_l = 0;
	this->n_macroparticles_lost = 0;
	this->id = new int[n_macroparticles];
	for (int i = 0; i < n_macroparticles; ++i) {
		id[i] = i + 1;
	}
}

inline int Beams::n_macroparticles_alive() {

	return n_macroparticles - n_macroparticles_lost;
}

inline void Beams::statistics() {
	ftype m_dE, m_dt, s_dE, s_dt;
	m_dt = m_dE = s_dE = s_dt = 0;
	int n = 0;
	for (int i = 0; i < n_macroparticles; ++i) {
		if (id[i] != 0) {
			m_dE += dE[i];
			m_dt += dt[i];
			n++;
		}
	}
	mean_dE = m_dE /= n;
	mean_dt = m_dt /= n;
	for (int i = 0; i < n_macroparticles; ++i) {
		if (id[i] != 0) {
			s_dE += (dE[i] - m_dE) * (dE[i] - m_dE);
			s_dt += (dt[i] - m_dt) * (dt[i] - m_dt);
		}
	}
	sigma_dE = sqrt(s_dE / n);
	sigma_dt = sqrt(s_dt / n);

	epsn_rms_l = pi * sigma_dE * sigma_dt; // in eVs

	//Losses
	n_macroparticles_lost = n_macroparticles - n;
}

inline void Beams::losses_longitudinal_cut(ftype dt_min, ftype dt_max) {
	for (int i = 0; i < n_macroparticles; ++i) {
		ftype a = (dt[i] - dt_min) * (dt_max - dt[i]);
		if (a < 0)
			id[i] = 0;
	}
}

inline void Beams::losses_energy_cut(ftype dE_min, ftype dE_max) {
	for (int i = 0; i < n_macroparticles; ++i) {
		ftype a = (dE[i] - dE_min) * (dE_max - dE[i]);
		if (a < 0)
			id[i] = 0;
	}
}

#endif /* BEAMS_BEAMS_H_ */
