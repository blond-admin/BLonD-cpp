/*
 * Monitors.cpp
 *
 *  Created on: August 22, 2016
 *      Author: kiliakis
 */

#include <blond/monitors/Monitors.h>
#include <string>
#include <iostream>

using namespace H5;

SlicesMonitor::SlicesMonitor(std::string filename, int n_turns, Slices *slices)
{
    fFileName = filename;
    fNTurns = n_turns;
    fSlices = slices;
    fITurn = 0;
    fFile = new H5File(fFileName, H5F_ACC_TRUNC);
    fGroup = new Group(fFile->createGroup("Slices"));
}

SlicesMonitor::~SlicesMonitor()
{
    delete fFile;
    delete fGroup;
}

void SlicesMonitor::track()
{
    if (fITurn == 0)
        create_data({fSlices->n_slices, (uint) fNTurns});
    write_data();
    fITurn++;
}


void SlicesMonitor::write_data()
{
    auto dataset = fGroup->openDataSet("n_macroparticles");
    hsize_t offset[2], count[2], stride[2], block[2];
    const uint size = fSlices->n_macroparticles.size();
    count[0] = size;
    count[1] = 1;
    offset[0] = 0;
    offset[1] = fITurn;
    stride[0] = stride[1] = 1;
    block[0] = block[1] = 1;
    DataSpace memspace(2, count, NULL);
    auto dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);

    auto temp = new int[size][1];
    for (uint i = 0; i < size; i++)
        temp[i][0] = fSlices->n_macroparticles[i];

    dataset.write(temp, PredType::NATIVE_INT, memspace, dataspace);
    delete[] temp;

}


void SlicesMonitor::create_data(const uint_vector_t dims)
{
    hsize_t dim[2], chunk[2];
    dim[0] = dims[0];
    dim[1] = dims[1];
    chunk[0] = dims[0];
    chunk[1] = 1;
    DataSpace dataspace(2, dim);
    DSetCreatPropList plist;
    plist.setChunk(2, chunk);
    plist.setDeflate(9);
    auto dset = DataSet(fGroup->createDataSet("n_macroparticles",
                        PredType::NATIVE_INT,
                        dataspace, plist));
}


void SlicesMonitor::close()
{
    fGroup->close();
    fFile->close();
}




BunchMonitor::BunchMonitor() {}
BunchMonitor::~BunchMonitor() {}
void BunchMonitor::track() {}
void BunchMonitor::init_data() {}
void BunchMonitor::init_buffer() {}
void BunchMonitor::write_buffer() {}
void BunchMonitor::write_data() {}
void BunchMonitor::open() {}
void BunchMonitor::close() {}
