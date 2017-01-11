
# Copyright 2016 CERN. This software is distributed under the
# terms of the GNU General Public Licence version 3 (GPL Version 3),
# copied verbatim in the file LICENCE.md.
# In applying this licence, CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.
# Project website: http://blond.web.cern.ch/

'''
**Module to plot simple stuff**

:Authors: **alasheen**
'''

from __future__ import division
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def plot_xy(xArray, yArray, figureName, dirName, saveFigName, xMax):

    plt.figure(figureName)
    plt.clf()
    plt.plot(xArray, yArray)
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    fign = dirName + '/'+saveFigName+'.png'
    if xMax >0:
        plt.xlim((0,xMax))
    plt.savefig(fign)



