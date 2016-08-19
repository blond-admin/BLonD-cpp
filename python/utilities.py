
# Copyright 2016 CERN. This software is distributed under the
# terms of the GNU General Public Licence version 3 (GPL Version 3),
# copied verbatim in the file LICENCE.md.
# In applying this licence, CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.
# Project website: http://blond.web.cern.ch/

'''

**Utilities to calculate Hamiltonian, separatrix, total voltage for the full ring.**

:Authors: **Danilo Quartullo**, **Helga Timko**, **Alexandre Lasheen**
'''


from __future__ import division
import warnings
import numpy as np


def separatrix(gp_n_sections, gp_charge, gp_t_rev,
               rfp_counter, rfp_voltage, rfp_omega_rf,
               rfp_phi_RF, rfp_eta_0, rfp_beta,
               rfp_energy, rfp_n_rf, rfp_harmonic_0,
               rfp_phi_S, rfp_E_increment, dt, total_voltage=None):
    """Single RF sinusoidal separatrix.
    For the time being, for single RF section only or from total voltage.
    Uses beta, energy averaged over the turn.
    To be generalized."""

    warnings.filterwarnings("once")

    if gp_n_sections > 1:
        warnings.warn(
            "WARNING: The separatrix is not yet properly computed for several sections!")

    # Import RF and ring parameters at this moment
    counter = rfp_counter
    voltage = gp_charge*rfp_voltage
    omega_RF = rfp_omega_rf
    phi_RF = rfp_phi_RF

    eta0 = rfp_eta_0
    beta_sq = rfp_beta**2
    energy = rfp_energy

    # Projects time array into the range [-T_RF/2+t_RF, T_RF/2+t_RF]
    # if below transition and into the range [t_RF, t_RF+T_RF] if above transition.
    # T_RF = 2*pi/omega_RF, t_RF = - phi_RF/omega_RF
    if eta0 < 0:
        dt = time_modulo(dt, (phi_RF[0] - np.pi)/omega_RF[0],
                         2.*np.pi/omega_RF[0])
    elif eta0 > 0:
        dt = time_modulo(dt, phi_RF[0]/omega_RF[0], 2.*np.pi/omega_RF[0])

    # Single-harmonic RF system
    if rfp_n_rf == 1:

        h0 = rfp_harmonic_0

        if total_voltage == None:
            V0 = voltage[0]
        else:
            V0 = total_voltage[counter]

        phi_s = rfp_phi_S
        phi_b = omega_RF[0]*dt + phi_RF[0]

        separatrix_array = np.sqrt(beta_sq*energy*V0/(np.pi*eta0*h0) *
                                   (-np.cos(phi_b) - np.cos(phi_s) +
                                    (np.pi - phi_s - phi_b)*np.sin(phi_s)))

    # Multi-harmonic RF system
    else:
        denergy = rfp_E_increment
        T0 = gp_t_rev
        index_voltage = np.min(np.where(voltage > 0)[0])
        T_RF0 = 2*np.pi/omega_RF[index_voltage]

        # Find unstable fixed point

        dt_ufp = np.linspace(-phi_RF[index_voltage]/omega_RF[index_voltage]
                             - T_RF0/1000, T_RF0
                             - phi_RF[index_voltage]/omega_RF[index_voltage]
                             + T_RF0/1000, 1002)

        if eta0 < 0:
            dt_ufp -= 0.5*T_RF0
        Vtot = np.zeros(len(dt_ufp))

        # Construct waveform
        for i in range(rfp_n_rf):
            temp = np.sin(omega_RF[i]*dt_ufp + phi_RF[i])
            Vtot += voltage[i]*temp
        Vtot -= denergy

        # Find zero crossings
        zero_crossings = np.where(np.diff(np.sign(Vtot)))[0]

        # Interpolate UFP
        if eta0 < 0:
            i = -1
            ind = zero_crossings[i]
            while (Vtot[ind+1] - Vtot[ind]) > 0:
                i -= 1
                ind = zero_crossings[i]
        else:
            i = 0
            ind = zero_crossings[i]
            while (Vtot[ind+1] - Vtot[ind]) < 0:
                i += 1
                ind = zero_crossings[i]
        dt_ufp = dt_ufp[ind] + Vtot[ind]/(Vtot[ind] - Vtot[ind+1]) * \
            (dt_ufp[ind+1] - dt_ufp[ind])

        # Construct separatrix
        Vtot = np.zeros(len(dt))
        for i in range(rfp_n_rf):
            Vtot += voltage[i]*(np.cos(omega_RF[i]*dt_ufp + phi_RF[i]) -
                                np.cos(omega_RF[i]*dt + phi_RF[i]))/omega_RF[i]

        separatrix_array = np.sqrt(2*beta_sq*energy/(eta0*T0) *
                                   (Vtot + denergy*(dt_ufp - dt)))

    return separatrix_array


def phase_modulo_above_transition(phi):
    '''
    *Projects a phase array into the range -Pi/2 to +3*Pi/2.*
    '''

    return phi - 2.*np.pi*np.floor(phi/(2.*np.pi))


def phase_modulo_below_transition(phi):
    '''
    *Projects a phase array into the range -Pi/2 to +3*Pi/2.*
    '''

    return phi - 2.*np.pi*(np.floor(phi/(2.*np.pi) + 0.5))


def time_modulo(dt, dt_offset, T):
    '''
    *Returns dt projected onto the desired interval.*
    '''

    return dt - T*np.floor((dt + dt_offset)/T)


def gaussian_filter1d(x, sigma, order, mode):
    from scipy import ndimage
    return ndimage.gaussian_filter1d(x, sigma=sigma, order=order, mode=mode)


def gradient(x, dist):
    return np.gradient(x, dist)
