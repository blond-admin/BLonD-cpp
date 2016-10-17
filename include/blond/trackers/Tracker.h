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
    // Imported fields
    RfParameters *rfp;
    Beams *beam;
    int &section_index;
    int &counter;
    double &length_ratio;
    double &section_length;
    f_vector_t &t_rev;
    int &n_rf;
    f_vector_t &beta;
    double &charge;
    f_vector_2d_t &harmonic;
    f_vector_2d_t &voltage;
    f_vector_2d_t &phi_noise;
    f_vector_2d_t &phi_rf;
    f_vector_t &phi_s;
    f_vector_2d_t &omega_rf;
    f_vector_t &eta_0;
    f_vector_t &eta_1;
    f_vector_t &eta_2;
    int_vector_t &sign_eta_0;
    int &alpha_order;

    // double elapsed_time;
    int_vector_t indices_right_outside;
    int_vector_t indices_inside_frame;
    int_vector_t indices_left_outside;

    f_vector_t acceleration_kick;
    f_vector_t fRfVoltage;
    solver_type solver;
    double dE_max;
    bool rf_kick_interp;
    bool periodicity;

    LHCNoiseFB *noiseFB;
    PhaseLoop *PL;
    Slices *slices;
    TotalInducedVoltage *totalInducedVoltage;
    f_vector_t fTotalVoltage;

    void set_periodicity();
    void kick(f_vector_t &beam_dt, f_vector_t &beam_dE, const int index);
    inline void kick(const double *__restrict beam_dt, double *__restrict beam_dE,
                     const int n_rf, const double *__restrict voltage,
                     const double *__restrict omega_RF,
                     const double *__restrict phi_RF, const int n_macroparticles,
                     const double acc_kick);
    void drift(f_vector_t &beam_dt, f_vector_t &beam_dE, const int index);
    inline void drift(double *__restrict beam_dt,
                      const double *__restrict beam_dE, const solver_type solver,
                      const double T0, const double length_ratio,
                      const int alpha_order, const double eta_zero,
                      const double eta_one, const double eta_two,
                      const double beta, const double energy,
                      const int n_macroparticles);

    void track();
    void rf_voltage_calculation(int turn, Slices *slices);

    inline void horizontal_cut();
    RingAndRfSection(RfParameters *RfP = Context::RfP,
                     Beams *Beam = Context::Beam,
                     solver_type solver = simple, PhaseLoop *PL = nullptr,
                     LHCNoiseFB *NoiseFB = nullptr, bool periodicity = false,
                     double dE_max = 0, bool rf_kick_interp = false,
                     Slices *Slices = nullptr,
                     TotalInducedVoltage *TotalInducedVoltage = nullptr)
        : section_index(RfP->section_index),
          counter(RfP->counter),
          length_ratio(RfP->length_ratio),
          section_length(RfP->section_length),
          t_rev(RfP->t_rev),
          n_rf(RfP->n_rf),
          beta(RfP->beta),
          charge(RfP->charge),
          harmonic(RfP->harmonic),
          voltage(RfP->voltage),
          phi_noise(RfP->phi_noise),
          phi_rf(RfP->phi_rf),
          phi_s(RfP->phi_s),
          omega_rf(RfP->omega_rf),
          eta_0(RfP->eta_0),
          eta_1(RfP->eta_1),
          eta_2(RfP->eta_2),
          sign_eta_0(RfP->sign_eta_0),
          alpha_order(RfP->alpha_order)
    {
        beam = Beam;
        rfp = RfP;
        // this->elapsed_time = 0;
        this->solver = solver;
        this->PL = PL;
        this->noiseFB = NoiseFB;
        this->periodicity = periodicity;
        this->dE_max = dE_max;
        this->rf_kick_interp = rf_kick_interp;
        this->slices = Slices;
        this->totalInducedVoltage = TotalInducedVoltage;

        this->acceleration_kick.resize(rfp->E_increment.size());
        for (uint i = 0; i < rfp->E_increment.size(); ++i)
            acceleration_kick[i] = -rfp->E_increment[i];

        if (solver != simple && solver != full) {
            std::cerr << "ERROR: Choice of longitudinal solver not recognized!\n"
                      << "Aborting...\n";
            exit(-1);
        }

        if (alpha_order > 1) solver = full;

        if (rf_kick_interp && Slices == NULL) {
            std::cerr << "ERROR: A slices object is needed to use the"
                      << " kick_interp option\n";
            exit(-1);
        }

    }
    ~RingAndRfSection() {};
};

class API FullRingAndRf {
private:
public:
    enum main_harmonic_t { lowest_freq = 0, highest_voltage = 1 };

    f_vector_t fPotentialWell;
    f_vector_t fPotentialWellCoordinates;
    f_vector_t fTotalVoltage;
    double fRingCircumference;
    double fRingRadius;
    std::vector<RingAndRfSection * > fRingList;

    void track();
    void potential_well_generation(const int turn = 0,
                                   const int n_points = 100000,
                                   const double option = lowest_freq,
                                   const double dt_margin_percent = 0.0);

    FullRingAndRf(std::vector<RingAndRfSection * > &RingList);
    ~FullRingAndRf();
};

#endif /* TRACKERS_TRACKER_H_ */
