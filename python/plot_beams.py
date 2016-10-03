
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
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from utilities import separatrix
import os

def plot_long_phase_space(rfp_counter, rfp_omega_RF_0, rfp_phi_RF_0,
                          beam_id, beam_dt, beam_dE, xmin, xmax, ymin,
                          ymax, xunit, sampling, separatrix_plot,
                          gp_n_sections, gp_charge, gp_t_rev,
                          rfp_voltage, rfp_omega_rf, rfp_phi_RF,
                          rfp_eta_0, rfp_beta, rfp_energy, rfp_n_rf,
                          rfp_harmonic_0, rfp_phi_S, rfp_E_increment,
                          histograms_plot, dirname, alpha,
                          ):
    """
    Plot of longitudinal phase space. Optional use of histograms and separatrix.
    Choice of units: xunit = s, rad.
    For large amount of data, use "sampling" to plot a fraction of the data.
    """
    # print beam_dE[:10]
    # print beam_dt[:10]
    # print beam_id[:10]
    # print "ok till here"
    # Conversion from particle arrival time to RF phase
    if xunit == 'rad':
        omega_RF = rfp_omega_RF_0
        # RFSectionParameters.omega_RF[0, RFSectionParameters.counter[0]]
        phi_RF = rfp_phi_RF_0
        # RFSectionParameters.phi_RF[0, RFSectionParameters.counter[0]]

    # Definitions for placing the axes
    left, width = 0.115, 0.63
    bottom, height = 0.115, 0.63
    bottom_h = left_h = left+width+0.03

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # Prepare plot
    fig = plt.figure(1)
    fig.set_size_inches(8, 8)
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # Main plot: longitudinal distribution
    indlost = np.where(beam_id[::sampling] == 0)[0]  # particles lost
    indalive = np.where(beam_id[::sampling] != 0)[0]  # particles transmitted
    if xunit == 's':
        axScatter.set_xlabel(r"$\Delta t$ [s]")
        axScatter.scatter(beam_dt[indalive], beam_dE[indalive],
                          s=1, edgecolor='none', alpha=alpha)
        axScatter.scatter(beam_dt[indlost], beam_dE[indlost], c='0.5',
                          s=1, edgecolor='none')
        
    elif xunit == 'rad':
        axScatter.set_xlabel(r"$\varphi$ [rad]")
        axScatter.scatter(omega_RF*beam_dt[indalive] + phi_RF,
                          beam_dE[indalive], s=1, edgecolor='none',
                          alpha=alpha)
        axScatter.scatter(omega_RF*beam_dt[indlost] + phi_RF,
                          beam_dE[indlost], c='0.5', s=1,
                          edgecolor='none')
    axScatter.set_ylabel(r"$\Delta$E [eV]")
    axScatter.yaxis.labelpad = 1

    axScatter.set_xlim(xmin, xmax)
    axScatter.set_ylim(ymin, ymax)

    axScatter.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.figtext(0.95, 0.95, '%d turns' % rfp_counter,
                fontsize=16, ha='right', va='center')

    # Separatrix

    if separatrix_plot:
        x_sep = np.linspace(xmin, xmax, 1000)
        if xunit == 's':
            y_sep = separatrix(gp_n_sections, gp_charge, gp_t_rev,
                               rfp_counter, rfp_voltage, rfp_omega_rf,
                               rfp_phi_RF, rfp_eta_0, rfp_beta,
                               rfp_energy, rfp_n_rf, rfp_harmonic_0,
                               rfp_phi_S, rfp_E_increment, x_sep)
        elif xunit == 'rad':
            y_sep = separatrix(gp_n_sections, gp_charge, gp_t_rev,
                               rfp_counter, rfp_voltage, rfp_omega_rf,
                               rfp_phi_RF, rfp_eta_0, rfp_beta,
                               rfp_energy, rfp_n_rf, rfp_harmonic_0,
                               rfp_phi_S, rfp_E_increment,
                               (x_sep - phi_RF)/omega_RF)
        axScatter.plot(x_sep, y_sep, 'r')
        axScatter.plot(x_sep, - y_sep, 'r')

    # Phase and momentum histograms
    if histograms_plot:
        xbin = (xmax - xmin)/200.
        xh = np.arange(xmin, xmax + xbin, xbin)
        ybin = (ymax - ymin)/200.
        yh = np.arange(ymin, ymax + ybin, ybin)
        # print beam_dt[::sampling]
        # print xh
        if xunit == 's':
            axHistx.hist(beam_dt[::sampling], bins=xh, histtype='step')
        elif xunit == 'rad':
            axHistx.hist(
                omega_RF*beam_dt[::sampling] + phi_RF, bins=xh, histtype='step')

        axHisty.hist(beam_dE[::sampling], bins=yh, histtype='step',
                     orientation='horizontal')

        axHistx.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        axHisty.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        axHistx.axes.get_xaxis().set_visible(False)
        axHisty.axes.get_yaxis().set_visible(False)
        axHistx.set_xlim(xmin, xmax)
        axHisty.set_ylim(ymin, ymax)
        labels = axHisty.get_xticklabels()
        for label in labels:
            label.set_rotation(-90)
    # print "hello from python!"

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    fign = dirname + '/long_distr_'"%d" % rfp_counter+'.png'
    plt.savefig(fign)
    plt.clf()


def plot_bunch_length_evol(time_step, sigma_dt, output_freq=1,
                           dirname='fig'):
    """
    Plot of r.m.s. 4-sigma bunch length [s] as a function of time.
    """

    # Time step of plotting
    # time_step = RFSectionParameters.counter[0]

    # Get bunch length data in metres or seconds
    if output_freq < 1:
        output_freq = 1
    ndata = int(time_step/output_freq)
    t = output_freq*np.arange(ndata)
    bl = np.array(sigma_dt[0:ndata], dtype=np.double)
    bl *= 4
    bl[time_step:] = np.nan

    # Plot
    fig = plt.figure(1)
    fig.set_size_inches(8, 6)
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(t, bl, '.')
    ax.set_xlabel(r"No. turns [T$_0$]")
    ax.set_ylabel(r"Bunch length, $\Delta t_{4\sigma}$ r.m.s. [s]")
    if time_step > 100000:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    fign = dirname + '/bunch_length.png'
    plt.savefig(fign)
    plt.clf()


def plot_bunch_length_evol_gaussian(time_step, bunch_length_gaussian,
                                    output_freq=1, dirname='fig'):
    """
    Plot of Gaussian 4-sigma bunch length [s] as a function of time.
    Requires slices.
    """

    # Time step of plotting
    # time_step = RFSectionParameters.counter[0]

    # Get bunch length data in metres or nanoseconds
    if output_freq < 1:
        output_freq = 1
    ndata = int(time_step/output_freq)
    t = output_freq*np.arange(ndata)

    bl = np.array(
        bunch_length_gaussian[0:ndata], dtype=np.double)

    bl[time_step:] = np.nan

    # Plot
    fig = plt.figure(1)
    fig.set_size_inches(8, 6)
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(t, bl, '.')
    ax.set_xlabel(r"No. turns [T$_0$]")
    ax.set_ylabel(r"Bunch length, $\Delta t_{4\sigma}$ Gaussian fit [s]")
    if time_step > 100000:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    fign = dirname + '/bunch_length_Gaussian.png'
    plt.savefig(fign)
    plt.clf()


def plot_position_evol(time_step, mean_dt, output_freq=1,
                       style='.', dirname='fig'):

    # Time step of plotting
    # time_step = RFSectionParameters.counter[0]

    # Get position data [s]
    if output_freq < 1:
        output_freq = 1
    ndata = int(time_step/output_freq)
    t = output_freq*np.arange(ndata)
    # print t
    pos = np.array(mean_dt[0:ndata], dtype=np.double)
    pos[time_step:] = np.nan
    # print pos
    # Plot
    fig = plt.figure(1)
    fig.set_size_inches(8, 6)
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(t, pos, style)
    ax.set_xlabel(r"No. turns [T$_0$]")
    ax.set_ylabel(r"Bunch mean position, $<\Delta t>$ [s]")
    if time_step > 100000:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    fign = dirname+'/bunch_mean_position.png'
    plt.savefig(fign)
    plt.clf()


def plot_energy_evol(time_step, mean_dE, output_freq=1, style='.',
                     dirname='fig'):

    # Time step of plotting
    # time_step = RFSectionParameters.counter[0]

    # Get position data in metres or nanoseconds
    if output_freq < 1:
        output_freq = 1
    ndata = int(time_step/output_freq)
    t = output_freq*np.arange(ndata)
    nrg = np.array(mean_dE[0:ndata], dtype=np.double)
    nrg[time_step:] = np.nan

    # Plot
    fig = plt.figure(1)
    fig.set_size_inches(8, 6)
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(t, nrg, style)
    ax.set_xlabel(r"No. turns [T$_0$]")
    ax.set_ylabel(r"Bunch mean energy, $<\Delta E>$ [eV]")
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    if time_step > 100000:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # Save plot
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    fign = dirname+'/bunch_mean_energy.png'
    plt.savefig(fign)
    plt.clf()


def plot_transmitted_particles(time_step, parts_alive, output_freq=1,
                               style='.', dirname='fig'):

    # Time step of plotting
    # time_step = RFSectionParameters.counter[0]

    # Get position data in metres or nanoseconds
    if output_freq < 1:
        output_freq = 1
    ndata = int(time_step/output_freq)
    t = output_freq*np.arange(ndata)
    prtcls = np.array(parts_alive[0:ndata],
                      dtype=np.double)
    prtcls[time_step:] = np.nan
    # print prtcls
    # Plot
    plt.figure(1, figsize=(8, 6))
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(t, prtcls, style)
    ax.set_xlabel(r"No. turns [T$_0$]")
    ax.set_ylabel(r"Transmitted macro-particles [1]")
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    if time_step > 100000:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    # Save plot
    fign = dirname+'/bunch_transmitted_particles.png'
    plt.savefig(fign)
    plt.clf()
