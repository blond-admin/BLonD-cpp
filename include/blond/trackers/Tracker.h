/*
 * Tracker.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#ifndef TRACKERS_TRACKER_H_
#define TRACKERS_TRACKER_H_

#include <blond/beams/Slices.h>
#include <blond/impedances/InducedVoltage.h>
#include <blond/llrf/LHCNoiseFB.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/utilities.h>


class API RingAndRfSection {

public:
    enum solver_type { simple, full };

    RfParameters *fRfP;

    ftype elapsed_time;
    int_vector_t indices_right_outside;
    int_vector_t indices_inside_frame;
    int_vector_t indices_left_outside;

    f_vector_t acceleration_kick;
    f_vector_t fRfVoltage;
    solver_type solver;
    ftype dE_max;
    bool rf_kick_interp;
    bool periodicity;

    LHCNoiseFB *noiseFB;
    PhaseLoop *PL;
    Slices *slices;
    TotalInducedVoltage *totalInducedVoltage;
    f_vector_t fTotalVoltage;

    void set_periodicity();
    void kick(f_vector_t &beam_dt, f_vector_t &beam_dE, const uint index);
    inline void kick(const ftype *__restrict beam_dt, ftype *__restrict beam_dE,
                     const int n_rf, const ftype *__restrict voltage,
                     const ftype *__restrict omega_RF,
                     const ftype *__restrict phi_RF, const int n_macroparticles,
                     const ftype acc_kick);
    void drift(f_vector_t &beam_dt, f_vector_t &beam_dE, const uint index);
    inline void drift(ftype *__restrict beam_dt,
                      const ftype *__restrict beam_dE, const solver_type solver,
                      const ftype T0, const ftype length_ratio,
                      const uint alpha_order, const ftype eta_zero,
                      const ftype eta_one, const ftype eta_two,
                      const ftype beta, const ftype energy,
                      const int n_macroparticles);

    void track();
    void rf_voltage_calculation(uint turn, Slices *slices);

    inline void horizontal_cut();
    RingAndRfSection(RfParameters *rfp = Context::RfP,
                     solver_type solver = simple, PhaseLoop *PL = NULL,
                     LHCNoiseFB *NoiseFB = NULL, bool periodicity = false,
                     ftype dE_max = 0, bool rf_kick_interp = false,
                     Slices *Slices = NULL,
                     TotalInducedVoltage *TotalInducedVoltage = NULL);
    ~RingAndRfSection();
};

class API FullRingAndRf {
private:
public:
    enum main_harmonic_t { lowest_freq = 0, highest_voltage = 1 };

    f_vector_t fPotentialWell;
    f_vector_t fPotentialWellCoordinates;
    f_vector_t fTotalVoltage;
    ftype fRingCircumference;
    ftype fRingRadius;
    std::vector<RingAndRfSection *> fRingList;

    void track();
    void potential_well_generation(const uint turn = 0,
                                   const uint n_points = 100000,
                                   const ftype option = lowest_freq,
                                   const ftype dt_margin_percent = 0.0);

    FullRingAndRf(std::vector<RingAndRfSection *> &RingList);
    ~FullRingAndRf();
};

#endif /* TRACKERS_TRACKER_H_ */
