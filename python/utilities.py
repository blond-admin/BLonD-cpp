
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


def potential_well_cut(theta_coord_array, potential_array):
    '''
    *Function to cut the potential well in order to take only the separatrix
    (several cases according to the number of min/max).*
    '''

    # Check for the min/max of the potential well
    minmax_positions, minmax_values = minmax_location(theta_coord_array,
                                                      potential_array)
    min_theta_positions = minmax_positions[0]
    max_theta_positions = minmax_positions[1]
    max_potential_values = minmax_values[1]
    n_minima = len(min_theta_positions)
    n_maxima = len(max_theta_positions)

    if n_minima == 0:
        raise RuntimeError('The potential well has no minima...')
    if n_minima > n_maxima and n_maxima == 1:
        raise RuntimeError(
            'The potential well has more minima than maxima, and only one maximum')
    if n_maxima == 0:
        print ('Warning: The maximum of the potential well could not be found... \
                You may reconsider the options to calculate the potential well \
                as the main harmonic is probably not the expected one. \
                You may also increase the percentage of margin to compute \
                the potentiel well. The full potential well will be taken')
    elif n_maxima == 1:
        if min_theta_positions[0] > max_theta_positions[0]:
            saved_indexes = (potential_array < max_potential_values[0]) * \
                            (theta_coord_array > max_theta_positions[0])
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
            if potential_array[-1] < potential_array[0]:
                raise RuntimeError('The potential well is not well defined. \
                                    You may reconsider the options to calculate \
                                    the potential well as the main harmonic is \
                                    probably not the expected one.')
        else:
            saved_indexes = (potential_array < max_potential_values[0]) * \
                            (theta_coord_array < max_theta_positions[0])
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
            if potential_array[-1] > potential_array[0]:
                raise RuntimeError('The potential well is not well defined. \
                                    You may reconsider the options to calculate \
                                    the potential well as the main harmonic is \
                                    probably not the expected one.')
    elif n_maxima == 2:
        lower_maximum_value = np.min(max_potential_values)
        higher_maximum_value = np.max(max_potential_values)
        lower_maximum_theta = max_theta_positions[
            max_potential_values == lower_maximum_value]
        higher_maximum_theta = max_theta_positions[
            max_potential_values == higher_maximum_value]
        if len(lower_maximum_theta) == 2:
            saved_indexes = (potential_array < lower_maximum_value) * \
                            (theta_coord_array > lower_maximum_theta[0]) * \
                            (theta_coord_array < lower_maximum_theta[1])
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
        elif min_theta_positions[0] > lower_maximum_theta:
            saved_indexes = (potential_array < lower_maximum_value) * \
                            (theta_coord_array > lower_maximum_theta) * \
                            (theta_coord_array < higher_maximum_theta)
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
        else:
            saved_indexes = (potential_array < lower_maximum_value) * \
                            (theta_coord_array < lower_maximum_theta) * \
                            (theta_coord_array > higher_maximum_theta)
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
    elif n_maxima > 2:
        left_max_theta = np.min(max_theta_positions)
        right_max_theta = np.max(max_theta_positions)
        left_max_value = max_potential_values[
            max_theta_positions == left_max_theta]
        right_max_value = max_potential_values[
            max_theta_positions == right_max_theta]
        separatrix_value = np.min([left_max_value, right_max_value])
        saved_indexes = (theta_coord_array > left_max_theta) * \
            (theta_coord_array < right_max_theta) * \
            (potential_array < separatrix_value)
        theta_coord_sep = theta_coord_array[saved_indexes]
        potential_well_sep = potential_array[saved_indexes]

    return theta_coord_sep, potential_well_sep


def minmax_location(x, f):
    '''
    *Function to locate the minima and maxima of the f(x) numerical function.*
    '''

    f_derivative = np.diff(f)
    x_derivative = x[0:-1] + (x[1]-x[0])/2
    f_derivative = np.interp(x, x_derivative, f_derivative)

    f_derivative_second = np.diff(f_derivative)
    f_derivative_second = np.interp(x, x_derivative, f_derivative_second)

    warnings.filterwarnings("ignore")
    f_derivative_zeros = np.unique(np.append(
        np.where(f_derivative == 0), np.where(f_derivative[1:]/f_derivative[0:-1] < 0)))

    min_x_position = (x[f_derivative_zeros[f_derivative_second[f_derivative_zeros] >
                                           0] + 1] + x[f_derivative_zeros[f_derivative_second[f_derivative_zeros] > 0]])/2
    max_x_position = (x[f_derivative_zeros[f_derivative_second[f_derivative_zeros] <
                                           0] + 1] + x[f_derivative_zeros[f_derivative_second[f_derivative_zeros] < 0]])/2

    min_values = np.interp(min_x_position, x, f)
    max_values = np.interp(max_x_position, x, f)

    warnings.filterwarnings("default")

    return [min_x_position, max_x_position], [min_values, max_values]


def gauss(x, *p):
    A, x0, sx = p
    return A*np.exp(-(x-x0)**2/2./sx**2)


def curve_fit(bin_centers, n_macroparticles, p0):
    from scipy.optimize import curve_fit
    # n_macroparticles = np.array(n_macroparticles, dtype=float)
    # print "bin_centers", bin_centers
    # print "n_macroparticles", n_macroparticles
    # print "p0", p0
    result = curve_fit(gauss, bin_centers, n_macroparticles, p0)[0]
    # print result[1]
    # print result[2]
    return result


def beam_profile_filter_chebyshev(n_macroparticles, resolution,
                                  n_slices, filter_option):
    '''
    *This routine is filtering the beam profile with a type II Chebyshev
    filter. The input is a library having the following structure and
    informations:*

    filter_option = {'type':'chebyshev', 'pass_frequency':pass_frequency,
                     'stop_frequency':stop_frequency,'gain_pass':gain_pass,
                     'gain_stop':gain_stop}

    *The function returns nCoefficients, the number of coefficients used
    in the filter. You can also add the following option to plot and return
    the filter transfer function:*

    filter_option = {..., 'transfer_function_plot':True}
    '''
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    # import os
    from scipy.signal import cheb2ord, cheby2, filtfilt, freqz
    from numpy.fft import fftfreq, fft

    noisyProfile = np.array(n_macroparticles, dtype=float)

    freqSampling = 1 / resolution
    nyqFreq = freqSampling / 2.

    frequencyPass = float(filter_option['pass_frequency']) / nyqFreq
    frequencyStop = float(filter_option['stop_frequency']) / nyqFreq
    gainPass = float(filter_option['gain_pass'])
    gainStop = float(filter_option['gain_stop'])
    # Compute the lowest order for a Chebyshev Type II digital filter
    nCoefficients, wn = cheb2ord(
        frequencyPass, frequencyStop, gainPass, gainStop)
    # Compute the coefficients a Chebyshev Type II digital filter
    b, a = cheby2(nCoefficients, gainStop, wn, btype='low')

    # NOTE 2 lines of code have been commented out!!
    # Apply the filter forward and backwards to cancel the group delay
    # macroparticles = filtfilt(b, a, noisyProfile)
    # macroparticles = np.ascontiguousarray(macroparticles)
    # print "n_macroparticles: ", macroparticles

    if 'transfer_function_plot' in filter_option and \
            filter_option['transfer_function_plot'].lower() == "true":
                # Plot the filter transfer function
        w, transferGain = freqz(b, a=a, worN=n_slices)
        transferFreq = w / np.pi * nyqFreq
        group_delay = - \
            np.diff(-np.unwrap(-np.angle(transferGain))) / - \
            np.diff(w*freqSampling)

        plt.figure()
        ax1 = plt.subplot(311)
        plt.plot(transferFreq, 20 * np.log10(abs(transferGain)))
        plt.ylabel('Magnitude [dB]')
        plt.subplot(312, sharex=ax1)
        plt.plot(transferFreq, np.unwrap(-np.angle(transferGain)))
        plt.ylabel('Phase [rad]')
        plt.subplot(313, sharex=ax1)
        plt.plot(transferFreq[:-1], group_delay)
        plt.ylabel('Group delay [s]')
        plt.xlabel('Frequency [Hz]')

        plt.savefig("filter_transfer_function.png")
        plt.clf()

        # Plot the bunch spectrum and the filter transfer function
        plt.figure()

        plt.plot(np.fft.fftfreq(n_slices, resolution),
                 20 * np.log10(np.abs(np.fft.fft(noisyProfile))))
        plt.xlabel('Frequency [Hz]')
        plt.twinx()
        plt.plot(transferFreq, 20 * np.log10(abs(transferGain)), 'r')
        plt.xlim(0, plt.xlim()[1])

        plt.savefig("bunch_spectrum_filter_tranfer_function.png")
        plt.clf()

        res = np.array([nCoefficients], dtype=float)
        res = np.append(res,
                        [transferFreq, transferGain.real, transferGain.imag])
        # print n_slices
        # print res.shape
        return res
    else:
        # print "I am about to return"
        res = np.array([nCoefficients], dtype=float)
        res = np.append(res, [b, a])
        # print res
        return res
