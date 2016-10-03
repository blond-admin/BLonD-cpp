
# Copyright 2016 CERN. This software is distributed under the
# terms of the GNU General Public Licence version 3 (GPL Version 3),
# copied verbatim in the file LICENCE.md.
# In applying this licence, CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.
# Project website: http://blond.web.cern.ch/

'''
**Module to plot different bunch features**

:Authors: **Helga Timko**, **Danilo Quartullo**
'''

from __future__ import division
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_impedance_vs_frequency(counter,
                                frequency_array,
                                total_impedance,
                                beam_spectrum,
                                beam_spectrum_freq,
                                freq_array_2d,
                                impedance_real_2d,
                                impedance_imag_2d,
                                option1="sum", option2="no_spectrum",
                                option3="freq_fft", style='-',
                                cut_left_right=None, cut_up_down=None,
                                dirname='fig'):
    """
    Plot of impedance vs frequency.
    """
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    frequency_array = np.array(frequency_array)
    total_impedance = np.array(total_impedance)
    freq_array_2d = np.array(freq_array_2d)
    impedance_real_2d = np.array(impedance_real_2d)
    impedance_imag_2d = np.array(impedance_imag_2d)

    if option1 == "sum":
        ax1 = plt.subplots()[1]
        ax1.plot(frequency_array, total_impedance.real, style)
        ax1.plot(frequency_array, total_impedance.imag, style)
        if option2 == "spectrum":
            ax2 = ax1.twinx()
            ax2.plot(beam_spectrum_freq, np.abs(beam_spectrum))
        fign = dirname + '/sum_imp_vs_freq_fft' "%d" % counter + '.png'
        # plt.show()
        plt.savefig(fign)
        plt.clf()

    elif option1 == "single":
        fig0 = plt.figure(0)
        ax0 = fig0.add_subplot(111)

        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(111)

        for x, real, imag in zip(freq_array_2d, impedance_real_2d,
                                 impedance_imag_2d):
            ax0.plot(x, real, style)
            ax0.set_xlim(cut_left_right)
            ax0.set_ylim(cut_up_down)
            ax1.plot(x, imag, style)
            ax1.set_xlim(cut_left_right)
            ax1.set_ylim(cut_up_down)

        fign1 = dirname + '/real_imp_vs_'+option3+'_' "%d" % counter + '.png'
        if option2 == "spectrum":
            ax2 = ax0.twinx()
            ax2.plot(beam_spectrum_freq, np.abs(beam_spectrum))
            ax2.set_xlim(cut_left_right)
            ax2.set_ylim(cut_up_down)
        ax0.set_xlabel("Frequency [Hz]")
        ax0.set_ylabel("Real impedance [Ohm]")

        plt.figure(0)
        plt.savefig(fign1)
        plt.clf()
        fign2 = dirname + '/imag_imp_vs_'+option3+'_' "%d" % counter + '.png'

        if option2 == "spectrum":
            ax3 = ax1.twinx()
            ax3.plot(beam_spectrum_freq, np.abs(beam_spectrum))
            ax3.set_xlim(cut_left_right)
            ax3.set_ylim(cut_up_down)

        plt.figure(1)
        plt.savefig(fign2)
        plt.clf()


def plot_induced_voltage_vs_bin_centers(counter, bin_centers, induced_voltage,
                                        style='-', dirname='fig'):
    """
    Plot of induced voltage vs bin centers.
    """

    fig0 = plt.figure(0)
    fig0.set_size_inches(8, 6)
    ax0 = plt.axes([0.15, 0.1, 0.8, 0.8])
    plt.plot(bin_centers, induced_voltage, style)
    ax0.set_xlabel("Time [s]")
    ax0.set_ylabel("Induced voltage [V]")

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    fign = dirname + '/induced_voltage_' "%d" % counter + '.png'
    plt.savefig(fign)
    plt.clf()
