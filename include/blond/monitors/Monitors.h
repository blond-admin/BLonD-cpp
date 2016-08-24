/*
 * monitors.h
 *
 *  Created on: August 22, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDE_BLOND_MONITORS_H_
#define INCLUDE_BLOND_MONITORS_H_

#include <H5Cpp.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <blond/configuration.h>
#include <blond/utilities.h>
#include <blond/beams/Beams.h> 
#include <blond/beams/Slices.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/llrf/LHCNoiseFB.h>


class API SlicesMonitor {
public:
    H5::H5File *fFile;
    H5::Group *fGroup;
    // H5::DataSet *fDataSet;
    std::string fFileName;

    Slices *fSlices;
    int fITurn;
    int fNTurns;

    void track();
    void write_data();
    void create_data(const uint_vector_t dims);
    void close();
    SlicesMonitor(std::string filename, int n_turns, Slices *slices);
    ~SlicesMonitor();
};

class API BunchMonitor {
public:
    H5::H5File *fH5File;
    H5::Group *fH5Group;

    std::string fFileName;
    int fNTurns;
    int fITurn;
    int fBufferTime;
    RfParameters *fRfP;
    Beams *fBeam;
    Slices *fSlices;
    PhaseLoop *fPL;
    LHCNoiseFB *fNoiseFB;
    bool fGaussian;
    void track();
    void init_data(const int dimension);
    void init_buffer();
    void write_buffer();
    void write_data();
    void open();
    void close();
    BunchMonitor(GeneralParameters *GP, RfParameters *RfP, Beams *Beam,
                 std::string filename, int buffer_time = 0,
                 Slices *Slices = NULL, PhaseLoop *PL = NULL,
                 LHCNoiseFB *noiseFB = NULL);
    ~BunchMonitor();
};

#endif /* INCLUDE_BLOND_MONITORS_H_ */
