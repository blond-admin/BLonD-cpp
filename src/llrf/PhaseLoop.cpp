/*
 * PhaseLoop.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: kiliakis
 */

#include "PhaseLoop.h"

PhaseLoop::PhaseLoop(ftype *PL_gain, ftype window_coefficient, int _delay,
                     ftype *_phaseNoise, ftype *_LHCNoiseFB) {

    // General initializations
    this->delay = _delay;
    this->alpha = window_coefficient;
    this->gain = PL_gain;
    this->RFnoise = _phaseNoise;
    this->noiseFB = _LHCNoiseFB;
    // End of general initializations
}

LHC::LHC(ftype *PL_gain, ftype SL_gain, ftype window_coefficient,
         ftype *_phaseNoise, ftype *_LHCNoiseFB, int _delay) {

    // PhaseLoop(PL_gain, window_coefficient, _delay, _phaseNoise, _LHCNoiseFB);
    // General Initializations
    this->delay = _delay;
    this->alpha = window_coefficient;
    this->gain = PL_gain;
    this->RFnoise = _phaseNoise;
    this->noiseFB = _LHCNoiseFB;
    // End of general initializations

    this->gain2 = SL_gain;
    this->lhc_y = 0;
    this->lhc_a = new ftype[GP->n_turns + 1];
    this->lhc_t = new ftype[GP->n_turns + 1];
    // this->domega_RF = new ftype[GP->n_turns + 1];
    // dump(gain, 10, "gain\n");
    if (gain2 != 0) {
        for (int i = 0; i < GP->n_turns + 1; ++i) {
            lhc_a[i] = 5.25 - RfP->omega_s0[i] / (constant::pi * 40);
        }
        for (int i = 0; i < GP->n_turns + 1; ++i) {
            lhc_t[i] = (2 * constant::pi * RfP->Qs[i] * sqrt(lhc_a[i])) /
                       sqrt(1 +
                            gain[i] / gain2 *
                            sqrt((1 + 1 / lhc_a[i]) / (1 + lhc_a[i])));
        }
    } else {
        util::zero(lhc_a, GP->n_turns + 1);
        util::zero(lhc_t, GP->n_turns + 1);
    }
}

LHC::~LHC() {
    util::delete_array(lhc_a);
    util::delete_array(lhc_t);
    util::delete_array(gain);
}

PSB::PSB(ftype *PL_gain, ftype *_RL_gain, ftype _PL_period, ftype _RL_period,
         ftype *_coefficients, ftype window_coefficient, ftype *_phaseNoise,
         ftype *_LHCNoiseFB, int _delay) {

    PhaseLoop(PL_gain, window_coefficient, _delay, _phaseNoise, _LHCNoiseFB);

    // figure out what to do with these pairs
    this->gain2 = new ftype[2];
    this->dt = new int[2];

    if (_RL_gain == NULL) {
        gain2[0] = 0;
        gain2[1] = 0;
    } else {
        gain2[0] = _RL_gain[0];
        gain2[1] = _RL_gain[1];
    }

    if (_PL_period == 0)
        dt[0] = 10e-6;
    else
        dt[0] = _PL_period;

    if (_RL_period == 0)
        dt[1] = 7;
    else
        dt[1] = _RL_period;

    this->PL_counter = 1;
    this->on_time.push_back(0);

    precalculate_time();

    // TODO coefficients
    // How many can i have?

    //*Memory of previous phase correction, for phase loop.*
    dphi_av = 0;
    dphi_av_prev = 0;
    //*Memory of previous relative radial correction, for rad loop.*
    drho_prev = 0;
    //*Accumulated time for radial loop*
    t_accum = 0;
    //*Phase loop frequency correction [1/s]*
    domega_PL = 0;
    //*Radial loop frequency correction [1/s]*
    domega_RL = 0;
    // domega_RF = 0;
}

PSB::~PSB() {
    on_time.clear();
    util::delete_array(coefficients);
    util::delete_array(gain);
    util::delete_array(gain2);
    util::delete_array(dt);
}

void PhaseLoop::beam_phase() {
    /*
     *Beam phase measured at the main RF frequency and phase. The beam is
     convolved with the window function of the band-pass filter of the
     machine. The coefficients of sine and cosine components determine the
     beam phase, projected to the range -Pi/2 to 3/2 Pi. Note that this beam
     phase is already w.r.t. the instantaneous RF phase.*
     */

    // Main RF frequency at the present turn
    // omega_RF = self.rf_params.omega_RF[0,self.rf_params.counter[0]]
    ftype omega_RF = RfP->omega_RF[RfP->counter];
    ftype phi_RF = RfP->phi_RF[RfP->counter];
    // Convolve with window function
    //
    ftype *base = new ftype[Slice->n_slices];
    ftype *array = new ftype[Slice->n_slices];

#pragma omp parallel for
    for (int i = 0; i < Slice->n_slices; ++i) {
        ftype a = alpha * Slice->bin_centers[i];
        base[i] = exp(a) * Slice->n_macroparticles[i];
    }

#pragma omp parallel for
    for (int i = 0; i < Slice->n_slices; ++i) {
        ftype a = omega_RF * Slice->bin_centers[i] + phi_RF;
        array[i] = base[i] * sin(a);
    }
    ftype scoeff =
            mymath::trapezoid(array, Slice->bin_centers, Slice->n_slices);

#pragma omp parallel for
    for (int i = 0; i < Slice->n_slices; ++i) {
        ftype a = omega_RF * Slice->bin_centers[i] + phi_RF;
        array[i] = base[i] * cos(a);
    }

    ftype ccoeff =
            mymath::trapezoid(array, Slice->bin_centers, Slice->n_slices);

    phi_beam = atan(scoeff / ccoeff) + constant::pi;

    delete[] base;
    delete[] array;
}

void PhaseLoop::phase_difference() {
    /*
     Phase difference between beam and RF phase of the main RF system.
     Optional: add RF phase noise through dphi directly.*
     */

    // Correct for design stable phase
    int counter = RfP->counter;
    dphi = phi_beam - RfP->phi_s[counter];
    // Possibility to add RF phase noise through the PL
    if (RFnoise != NULL) {
        if (noiseFB != NULL) {
        }
    }
}

void PhaseLoop::radial_steering_from_freq() { }

void PhaseLoop::default_track() {
    /*
     Calculate PL correction on main RF frequency depending on machine.
     Update the RF phase and frequency of the next turn for all systems.
     */

    int counter = RfP->counter + 1;
    int turns = GP->n_turns;
    // Update the RF frequency of all systems for the next turn
    for (int i = 0; i < RfP->n_rf; ++i) {
        int row = i * (turns + 1);
        RfP->omega_RF[row + counter] +=
                domega_RF * RfP->harmonic[row + counter] / RfP->harmonic[counter];
    }

    // Update the RF phase of all systems for the next turn
    // Accumulated phase offset due to PL in each RF system
    for (int i = 0; i < RfP->n_rf; ++i) {
        int row = i * (turns + 1);
        RfP->dphi_RF[i] +=
                2 * constant::pi * RfP->harmonic[row + counter] *
                (RfP->omega_RF[row + counter] - RfP->omega_RF_d[row + counter]) /
                RfP->omega_RF_d[row + counter];
    }

    // Total phase offset
    for (int i = 0; i < RfP->n_rf; ++i) {
        int row = i * (turns + 1);
        RfP->phi_RF[row + counter] += RfP->dphi_RF[i];
    }
}

void LHC::track() {
    /*
     Calculation of the LHC RF frequency correction from the phase difference
     between beam and RF (actual synchronous phase). The transfer function is

     .. math::
     \\Delta \\omega_{rf}^{PL} = - g_{PL} (\\Delta\\varphi_{PL} + \\phi_{N})

     where the phase noise for the controlled blow-up can be optionally
     activated.
     Using 'gain2', a synchro loop can be activated in addition to remove
     long-term frequency drifts:

     .. math::
     \\Delta \\omega_{rf}^{SL} = - g_{SL} (y + a \\Delta\\varphi_{rf}) ,

     where we use the recursion

     .. math::
     y_{n+1} = (1 - \\tau) y_n + (1 - a) \\tau \\Delta\\varphi_{rf} ,

     with a and \tau being defined through the synchrotron frequency f_s and
     the synchrotron tune Q_s as

     .. math::
     a (f_s) \\equiv 5.25 - \\frac{f_s}{\\pi 40~\\text{Hz}} ,

     .. math::
     \\tau(f_s) \\equiv 2 \\pi Q_s \\sqrt{ \\frac{a}{1 + \\frac{g_{PL}}{g_{SL}}
     \\sqrt{\\frac{1 + 1/a}{1 + a}} }}
     */
    int counter = RfP->counter;
    ftype dphi_RF = RfP->dphi_RF[0];

    beam_phase();
    phase_difference();

    // Frequency correction from phase loop and synchro loop

    domega_RF = -gain[counter] * dphi -
                gain2 * (lhc_y + lhc_a[counter] * (dphi_RF + reference));

    // Update recursion variable
    lhc_y = (1 - lhc_t[counter]) * lhc_y +
            (1 - lhc_a[counter]) * lhc_t[counter] * (dphi_RF + reference);

    default_track();
}

void PSB::track() {
    /*
     Calculation of the PSB RF frequency correction from the phase difference
     between beam and RF (actual synchronous phase). The transfer function is

     .. math::
     \\Delta \\omega_{RF} = g(t) \\frac{a_0 \\Delta\\Phi_{PL}^2 + a_1
     \\Delta\\Phi_{PL} + a_2 }{b_0 \\Delta\\Phi_{PL}^2 + b_1 \\Delta\\Phi_{PL} +
     b_2}

     Input g through gain and [a_0, a_1, a_2, b_0, b_1, b_2] through
     coefficients.
     */

    // Average phase error while frequency is updated
    int counter = RfP->counter;

    beam_phase();
    phase_difference();
    dphi_av += dphi;
    t_accum += GP->t_rev[counter];

    // dprintf("Before if counter is %d\n", counter);
    // Phase loop active on certain turns

    if (counter == on_time[PL_counter] && counter > delay) {
        dphi_av /= (on_time[PL_counter] - on_time[PL_counter - 1]);
        domega_PL = 0.998 * domega_PL -
                    gain[counter] * (dphi_av - dphi_av_prev + reference);
    }
    // dprintf("After if\n");

    // Update averaging variables
    dphi_av_prev = dphi_av;
    dphi_av = 0;
    // Add correction from radial loop
    if (PL_counter % dt[1] == 0) {
        drho = (RfP->omega_RF[counter] - RfP->omega_RF_d[counter]) /
               (RfP->omega_RF_d[counter] *
                (1 / (GP->alpha[0] * RfP->gamma(counter) *
                      RfP->gamma(counter)) -
                 1)) +
               reference;

        domega_RL = domega_RL - gain2[0] * (drho - drho_prev) -
                    gain2[1] * drho * t_accum;

        drho_prev = drho;
        t_accum = 0;
    }

    // Counter to pick the next time step when the PL will be active
    PL_counter++;
    // Apply frequency correction
    domega_RF = domega_PL + domega_RL;
}

// TODO test this function
void PSB::radial_difference() {
    // Radial difference between beam and design orbit.*
    int counter = RfP->counter;
    int n = 0;
    ftype sum = 0;
    // ftype array[GP->n_turns];
    for (int i = 0; i < GP->n_turns; ++i) {
        if (Beam->dt[i] > Slice->bin_centers[0] &&
            Beam->dt[i] < Slice->bin_centers[Slice->n_slices - 1]) {
            sum += Beam->dE[i];
            n++;
        }
    }
    average_dE = sum / n;
    drho = GP->alpha[0] * GP->ring_radius * average_dE /
           (GP->beta[counter] * GP->beta[counter] * GP->energy[counter]);
}

// TODO Test this function
void PSB::precalculate_time() {
    /*
     For machines like the PSB, where the PL acts only in certain time
     intervals, pre-calculate on which turns to act.
     */
    unsigned int n = delay + 1;

    while (n < GP->t_rev.size()) {
        dprintf("dt[0] = %d\n", dt[0]);
        ftype summa = 0;
        while (summa < dt[0]) {
            try {
                summa += GP->t_rev[n];
                n++;
            } catch (...) {
                on_time.push_back(0);
                return;
            }
        }
        on_time.push_back(n - 1);
    }
}
