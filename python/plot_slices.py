
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


def plot_beam_profile(bin_centers, n_macroparticles, counter, style, dirname):
    """
    Plot of longitudinal beam profile
    """

    fig = plt.figure(1)
    fig.set_size_inches(8, 6)
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(bin_centers, n_macroparticles, style)

    ax.set_xlabel(r"$\Delta t$ [s]")
    ax.set_ylabel('Beam profile [arb. units]')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.figtext(0.95, 0.95, '%d turns' % counter, fontsize=16, ha='right',
                va='center')

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    fign = dirname + '/beam_profile_' "%d" % counter + '.png'
    plt.savefig(fign)
    plt.clf()


def plot_beam_profile_derivative(x, derivative, counter, style,
                                 dirname, modes):
    """
    Plot of the derivative of the longitudinal beam profile.
    """
    modes = modes.split(' ')
    for a, b, m in zip(x, derivative, modes):
        plt.plot(a, b, style, label=m)

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    fign = dirname + '/beam_profile_derivative_' "%d" % counter + '.png'
    plt.legend()
    plt.savefig(fign)
    plt.clf()


def plot_beam_spectrum(beam_spectrum_freq, beam_spectrum,
                       counter, style, dirname):
    """
    Plot of longitudinal beam profile
    """

    plt.figure(1, figsize=(8, 6))
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(beam_spectrum_freq, np.absolute(
        beam_spectrum), style)

    ax.set_xlabel(r"Frequency [Hz]")
    ax.set_ylabel('Beam spectrum, absolute value [arb. units]')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_xlim(0, 5.e9)

    plt.figtext(0.95, 0.95, '%d turns' % counter, fontsize=16, ha='right',
                va='center')

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    fign = dirname + '/beam_spectrum_' "%d" % counter + '.png'
    plt.savefig(fign)
    plt.clf()
