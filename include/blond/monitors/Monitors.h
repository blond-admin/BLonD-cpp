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


void *read_1D(std::string fname, std::string dsname,
              std::string type, hsize_t dims[]);

void *read_2D(std::string fname, std::string dsname,
              std::string type, hsize_t dims[]);


class API SlicesMonitor {
public:
    H5::H5File *fFile;
    H5::Group *fGroup;
    // H5::DataSet *fDataSet;
    std::string fFileName;

    Slices *fSlices;
    int fITurn;
    int fNTurns;
    int fCompressionLevel;

    void track();
    void write_data();
    void create_data(const uint_vector_t dims);
    void close();
    SlicesMonitor(std::string filename, int n_turns, Slices *slices,
                  int compression_level = 9);
    ~SlicesMonitor();
};

class API BunchMonitor {
private:
    int_vector_t b_np_alive;
    f_vector_t b_mean_dt;
    f_vector_t b_mean_dE;
    f_vector_t b_sigma_dt;
    f_vector_t b_sigma_dE;
    f_vector_t b_epsn_rms;
    f_vector_t b_bl_gauss;
    f_vector_t b_PL_omegaRF;
    f_vector_t b_PL_phiRF;
    f_vector_t b_PL_bunch_phase;
    f_vector_t b_PL_phase_corr;
    f_vector_t b_PL_omegaRF_corr;
    f_vector_t b_SL_dphiRF;
    f_vector_t b_RL_drho;
    f_vector_t b_LHCnoiseFB_factor;
    f_vector_t b_LHCnoiseFB_bl;
    f_vector_2d_t b_LHCnoiseFB_bl_bbb;
public:
    H5::H5File *fH5File;
    H5::Group *fH5Group;

    std::string fFileName;
    int fNTurns;
    int fITurn;
    int fBufferTime;
    int fCompressionLevel;
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
    // void open();
    void close();
    BunchMonitor(GeneralParameters *GP, RfParameters *RfP, Beams *Beam,
                 std::string filename, int buffer_time = 0,
                 Slices *Slices = NULL, PhaseLoop *PL = NULL,
                 LHCNoiseFB *noiseFB = NULL, int compression_level = 9);
    ~BunchMonitor();
};

#endif /* INCLUDE_BLOND_MONITORS_H_ */
