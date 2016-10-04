/*
 * Beams.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include <blond/beams/Beams.h>
#include <blond/constants.h>
#include <blond/math_functions.h>
#include <blond/trackers/utilities.h>

Beams::Beams(const uint n_macroparticles, const long long intensity)
{

    this->n_macroparticles = n_macroparticles;
    this->intensity = intensity;
    this->dt.resize(n_macroparticles);
    this->dE.resize(n_macroparticles);
    this->id.resize(n_macroparticles, 1);
    // this->id = mymath::arange<int>(1, n_macroparticles + 1);
    this->mean_dt = this->mean_dE = 0;
    this->sigma_dt = this->sigma_dE = 0;
    this->ratio = intensity / n_macroparticles;
    this->epsn_rms_l = 0;
    this->n_macroparticles_lost = 0;
}

Beams::~Beams() {}

uint Beams::n_macroparticles_alive()
{
    return n_macroparticles - n_macroparticles_lost;
}

void Beams::statistics()
{
    statistics(dE.data(), dt.data(), id.data(), dE.size());
}

void Beams::statistics(const double *__restrict dE,
                       const double *__restrict dt,
                       const int *__restrict id,
                       const int size)
{
    double m_dE, m_dt, s_dE, s_dt;
    m_dt = m_dE = s_dE = s_dt = 0.0;
    int n = 0;

    #pragma omp parallel for reduction(+:m_dE, m_dt, n)
    for (int i = 0; i < size; ++i) {
        m_dE += id[i] * dE[i];
        m_dt += id[i] * dt[i];
        n += id[i];
    }

    mean_dE = m_dE /= n;
    mean_dt = m_dt /= n;

    #pragma omp parallel for reduction(+:s_dE, s_dt)
    for (int i = 0; i < size; ++i) {
        s_dE += id[i] * (dE[i] - mean_dE) * (dE[i] - mean_dE);
        s_dt += id[i] * (dt[i] - mean_dt) * (dt[i] - mean_dt);
    }
    sigma_dE = std::sqrt(s_dE / n);
    sigma_dt = std::sqrt(s_dt / n);
    epsn_rms_l = constant::pi * sigma_dE * sigma_dt; // in eVs
    // Losses
    n_macroparticles_lost = n_macroparticles - n;

}


void Beams::losses_longitudinal_cut(const ftype dt_min, const ftype dt_max)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)n_macroparticles; i++)
        id[i] = (dt[i] - dt_min) * (dt_max - dt[i]) < 0 ? 0 : id[i];
}

void Beams::losses_energy_cut(const ftype dE_min, const ftype dE_max)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)n_macroparticles; ++i)
        id[i] = (dE[i] - dE_min) * (dE_max - dE[i]) < 0 ? 0 : id[i];
}


void Beams::losses_separatrix(GeneralParameters *GP, RfParameters *RfP)
{
    auto index = is_in_separatrix(GP, RfP, this, dt, dE);
    #pragma omp parallel for
    for (int i = 0; i < (int) n_macroparticles; i++)
        id[i] = id[i] * index[i];
}
/*
void Beams::losses_longitudinal_cut(const ftype* __restrict dt,
                                    const ftype dt_min, const ftype dt_max,
                                    int* __restrict id) {

#pragma omp parallel for
    for (int i = 0; i < (int)n_macroparticles; i++) {
        id[i] = (dt[i] - dt_min) * (dt_max - dt[i]) < 0 ? 0 : id[i];
    }
}

void Beams::losses_energy_cut(const ftype* __restrict dE, const ftype dE_min,
                              const ftype dE_max, int* __restrict id) {
#pragma omp parallel for
    for (int i = 0; i < (int)n_macroparticles; ++i) {
        id[i] = (dE[i] - dE_min) * (dE_max - dE[i]) < 0 ? 0 : id[i];
    }
}
*/
