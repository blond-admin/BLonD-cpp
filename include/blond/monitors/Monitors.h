/*
 * monitors.h
 *
 *  Created on: August 22, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDE_BLOND_MONITORS_H_
#define INCLUDE_BLOND_MONITORS_H_

#include <blond/configuration.h>
#include <blond/utilities.h>
#include <blond/beams/Slices.h>
#include <H5Cpp.h>


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
    void track();
    void init_data();
    void init_buffer();
    void write_buffer();
    void write_data();
    void open();
    void close();
    BunchMonitor();
    ~BunchMonitor();
};

#endif /* INCLUDE_BLOND_MONITORS_H_ */
