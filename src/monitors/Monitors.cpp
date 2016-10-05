/*
 * Monitors.cpp
 *
 *  Created on: August 22, 2016
 *      Author: kiliakis
 */

#include <blond/monitors/Monitors.h>
#include <string>
#include <iostream>
#include <algorithm>
using namespace H5;


void *read_2D(std::string fname, std::string dsname,
              std::string type, hsize_t dims[])
{
    auto file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5T_class_t class_id;
    size_t type_size;
    H5LTget_dataset_info(file_id, dsname.c_str(),
                         dims, &class_id, &type_size);
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    void *res = NULL;

    if (type == "int" && sizeof(int) == type_size) {
        res = malloc(dims[0] * dims[1] * type_size);
        H5LTread_dataset_int(file_id, dsname.c_str(), (int *)res);
    } else if (type == "double" && sizeof(double) == type_size) {
        res = malloc(dims[0] * dims[1] * type_size);
        H5LTread_dataset_double(file_id, dsname.c_str(), (double *)res);
    } else if (type == "float" && sizeof(float) == type_size) {
        res = malloc(dims[0] * dims[1] * type_size);
        H5LTread_dataset_float(file_id, dsname.c_str(), (float *)res);
    } else {
        std::cerr << "[Monitors] Requested type not supported or not compatible\n"
                  << "with actual data type size\n";
    }
    return res;
}

void *read_1D(std::string fname, std::string dsname,
              std::string type, hsize_t dims[])
{
    auto file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5T_class_t class_id;
    size_t type_size;
    H5LTget_dataset_info(file_id, dsname.c_str(),
                         dims, &class_id, &type_size);
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    void *res = NULL;

    if (type == "int" && sizeof(int) == type_size) {
        res = malloc(dims[0] * type_size);
        H5LTread_dataset_int(file_id, dsname.c_str(), (int *)res);
    } else if (type == "double" && sizeof(double) == type_size) {
        res = malloc(dims[0] * type_size);
        H5LTread_dataset_double(file_id, dsname.c_str(), (double *)res);
    } else if (type == "float" && sizeof(float) == type_size) {
        res = malloc(dims[0] * type_size);
        H5LTread_dataset_float(file_id, dsname.c_str(), (float *)res);
    } else {
        std::cerr << "[Monitors] Requested type not supported or not compatible\n"
                  << "with actual data type size\n";
    }
    return res;
}

SlicesMonitor::SlicesMonitor(std::string filename, int n_turns, Slices *slices,
                             int compression_level)
{
    fFileName = filename;
    fNTurns = n_turns;
    fSlices = slices;
    fITurn = 0;
    fFile = new H5File(fFileName, H5F_ACC_TRUNC);
    fGroup = new Group(fFile->createGroup("Slices"));
    fCompressionLevel = compression_level;
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
    plist.setDeflate(fCompressionLevel);
    auto dset = DataSet(fGroup->createDataSet("n_macroparticles",
                        PredType::NATIVE_INT,
                        dataspace, plist));
}


void SlicesMonitor::close()
{
    fGroup->close();
    fFile->close();
}




BunchMonitor::BunchMonitor(GeneralParameters *GP, RfParameters *RfP, Beams *Beam,
                           std::string filename, int buffer_time,
                           Slices *Slices, PhaseLoop *PL, LHCNoiseFB *noiseFB,
                           int compression_level)
{
    fFileName = filename;
    fNTurns = GP->n_turns;
    fITurn = 0;
    fBufferTime = buffer_time == 0 ? fNTurns : buffer_time;
    fRfP = RfP;
    fBeam = Beam;
    fSlices = Slices;
    fNoiseFB = noiseFB;
    fPL = PL;
    fH5File = new H5File(fFileName, H5F_ACC_TRUNC);
    fH5Group = new Group(fH5File->createGroup("Beam"));
    fCompressionLevel = compression_level;


    if (fSlices != NULL && fSlices->fit_option == fit_type::gaussian_fit)
        fGaussian = true;
    else
        fGaussian = false;

    init_data(fNTurns + 1);
    init_buffer();
    track();
}


BunchMonitor::~BunchMonitor()
{
    close();
}


void BunchMonitor::track()
{
    fBeam->statistics();
    if (fITurn <= fNTurns)
        write_buffer();
    fITurn++;
    if (fITurn > 0 && fITurn % fBufferTime == 0) {
        write_data();
        // init_buffer();
    }
}


void BunchMonitor::init_data(int dimension)
{
    int_vector_t intV(dimension);
    std::vector<real_t> realV(dimension);
    std::vector<double> doubleV(dimension);

    hsize_t dim[1], chunk[1];
    dim[0] = dimension;
    chunk[0] = dimension / 1000;

    DataSpace dataspace(1, dim);
    DSetCreatPropList plist;
    plist.setDeflate(fCompressionLevel);
    plist.setChunk(1, chunk);

    fBeam->statistics();


    auto dataset = fH5Group->createDataSet("n_macroparticles_alive",
                                           PredType::NATIVE_INT,
                                           dataspace, plist);
    intV[0] = fBeam->n_macroparticles_alive();
    dataset.write(intV.data(), PredType::NATIVE_INT);


    dataset = fH5Group->createDataSet("mean_dt",
                                      NATIVE_REAL_T,
                                      dataspace, plist);
    realV[0] = fBeam->mean_dt;
    dataset.write(realV.data(), NATIVE_REAL_T);


    dataset = fH5Group->createDataSet("mean_dE",
                                      NATIVE_REAL_T,
                                      dataspace, plist);
    realV[0] = fBeam->mean_dE;
    dataset.write(realV.data(), NATIVE_REAL_T);



    dataset = fH5Group->createDataSet("sigma_dt",
                                      NATIVE_REAL_T,
                                      dataspace, plist);
    realV[0] = fBeam->sigma_dt;
    dataset.write(realV.data(), NATIVE_REAL_T);


    dataset = fH5Group->createDataSet("sigma_dE",
                                      NATIVE_REAL_T,
                                      dataspace, plist);
    realV[0] = fBeam->sigma_dE;
    dataset.write(realV.data(), NATIVE_REAL_T);


    dataset = fH5Group->createDataSet("epsn_rms_l",
                                      NATIVE_REAL_T,
                                      dataspace, plist);
    realV[0] = fBeam->epsn_rms_l;
    dataset.write(realV.data(), NATIVE_REAL_T);


    if (fGaussian) {
        dataset = fH5Group->createDataSet("bunch_length_gaussian",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fSlices->bl_gauss;
        dataset.write(realV.data(), NATIVE_REAL_T);
    }


    if (fPL != NULL) {
        dataset = fH5Group->createDataSet("PL_omegaRF",
                                          PredType::NATIVE_DOUBLE,
                                          dataspace, plist);
        doubleV[0] = fRfP->omega_RF[0][0];
        dataset.write(doubleV.data(), PredType::NATIVE_DOUBLE);


        dataset = fH5Group->createDataSet("PL_phiRF",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fRfP->phi_RF[0][0];
        dataset.write(realV.data(), NATIVE_REAL_T);


        dataset = fH5Group->createDataSet("PL_bunch_phase",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fPL->phi_beam;
        dataset.write(realV.data(), NATIVE_REAL_T);


        dataset = fH5Group->createDataSet("PL_phase_corr",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fPL->dphi;
        dataset.write(realV.data(), NATIVE_REAL_T);


        dataset = fH5Group->createDataSet("PL_omegaRF_corr",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fPL->domega_RF;
        dataset.write(realV.data(), NATIVE_REAL_T);


        dataset = fH5Group->createDataSet("SL_dphiRF",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fRfP->dphi_RF[0];
        dataset.write(realV.data(), NATIVE_REAL_T);


        dataset = fH5Group->createDataSet("RL_drho",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fPL->drho;
        dataset.write(realV.data(), NATIVE_REAL_T);

    }


    if (fNoiseFB != NULL) {
        dataset = fH5Group->createDataSet("LHC_noise_FB_factor",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fNoiseFB->fX;
        dataset.write(realV.data(), NATIVE_REAL_T);


        dataset = fH5Group->createDataSet("LHC_noise_FB_bl",
                                          NATIVE_REAL_T,
                                          dataspace, plist);
        realV[0] = fNoiseFB->fBlMeas;
        dataset.write(realV.data(), NATIVE_REAL_T);

        if (!fNoiseFB->fBlMeasBBB.empty()) {
            hsize_t dim[2], chunk[2];
            dim[1] = fNoiseFB->fBlMeasBBB.size();
            dim[0] = dimension;
            chunk[0] = dimension / 1000;
            chunk[1] = 1;

            DataSpace dataspace(2, dim);
            DSetCreatPropList plist;
            plist.setChunk(2, chunk);
            plist.setDeflate(fCompressionLevel);

            dataset = fH5Group->createDataSet("LHC_noise_FB_bl_bbb",
                                              NATIVE_REAL_T,
                                              dataspace, plist);

            // realV = fNoiseFB->fBlMeasBBB;
            realV = std::vector<real_t>(fNoiseFB->fBlMeasBBB.begin(),
                                        fNoiseFB->fBlMeasBBB.end());
            dataset.write(realV.data(), NATIVE_REAL_T);
        }
    }

}

// TODO error when chosing buffer_time == 1,2
void BunchMonitor::write_data()
{
    // const int rank = dims.size();
    hsize_t offset[1], count[1], stride[1], block[1];
    count[0] = fBufferTime;
    offset[0] = fITurn - fBufferTime;
    stride[0] = 1;
    block[0] = 1;
    DataSpace memspace(1, count, NULL);
    // std::cout << fITurn << "\n";
    auto dataset = fH5Group->openDataSet("n_macroparticles_alive");
    auto dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dataset.write(b_np_alive.data(), PredType::NATIVE_INT,
                  memspace, dataspace);


    dataset = fH5Group->openDataSet("mean_dt");
    dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dataset.write(b_mean_dt.data(), NATIVE_REAL_T,
                  memspace, dataspace);


    dataset = fH5Group->openDataSet("mean_dE");
    dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dataset.write(b_mean_dE.data(), NATIVE_REAL_T,
                  memspace, dataspace);


    dataset = fH5Group->openDataSet("sigma_dt");
    dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dataset.write(b_sigma_dt.data(), NATIVE_REAL_T,
                  memspace, dataspace);


    dataset = fH5Group->openDataSet("sigma_dE");
    dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dataset.write(b_sigma_dE.data(), NATIVE_REAL_T,
                  memspace, dataspace);


    dataset = fH5Group->openDataSet("epsn_rms_l");
    dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dataset.write(b_epsn_rms.data(), NATIVE_REAL_T,
                  memspace, dataspace);

    if (fGaussian) {
        dataset = fH5Group->openDataSet("bunch_length_gaussian");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_bl_gauss.data(), NATIVE_REAL_T,
                      memspace, dataspace);

    }

    if (fPL != NULL) {
        dataset = fH5Group->openDataSet("PL_omegaRF");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_PL_omegaRF.data(), PredType::NATIVE_DOUBLE,
                      memspace, dataspace);

        dataset = fH5Group->openDataSet("PL_phiRF");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_PL_phiRF.data(), NATIVE_REAL_T,
                      memspace, dataspace);

        dataset = fH5Group->openDataSet("PL_bunch_phase");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_PL_bunch_phase.data(), NATIVE_REAL_T,
                      memspace, dataspace);

        dataset = fH5Group->openDataSet("PL_phase_corr");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_PL_phase_corr.data(), NATIVE_REAL_T,
                      memspace, dataspace);

        dataset = fH5Group->openDataSet("PL_omegaRF_corr");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_PL_omegaRF_corr.data(), NATIVE_REAL_T,
                      memspace, dataspace);


        dataset = fH5Group->openDataSet("SL_dphiRF");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_SL_dphiRF.data(), NATIVE_REAL_T,
                      memspace, dataspace);

        dataset = fH5Group->openDataSet("RL_drho");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_RL_drho.data(), NATIVE_REAL_T,
                      memspace, dataspace);
    }

    if (fNoiseFB != NULL) {
        dataset = fH5Group->openDataSet("LHC_noise_FB_factor");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_LHCnoiseFB_factor.data(), NATIVE_REAL_T,
                      memspace, dataspace);


        dataset = fH5Group->openDataSet("LHC_noise_FB_bl");
        dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        dataset.write(b_LHCnoiseFB_bl.data(), NATIVE_REAL_T,
                      memspace, dataspace);


        // TODO: Particularly test this part
        if (!fNoiseFB->fBlMeasBBB.empty()) {

            hsize_t offset[2], count[2], stride[2], block[2];
            const uint size = fNoiseFB->fBlMeasBBB.size();
            count[0] = 1;
            count[1] = size;
            offset[0] = fITurn;
            offset[1] = 0;
            stride[0] = stride[1] = 1;
            block[0] = block[1] = 1;

            DataSpace memspace(2, count, NULL);

            real_t *temp = new real_t[fBufferTime * size];
            for (int i = 0; i < fBufferTime; i++)
                for (uint j = 0; j < size; j++)
                    temp[i * size + j] = b_LHCnoiseFB_bl_bbb[i][j];


            dataset = fH5Group->openDataSet("LHC_noise_FB_bl_bbb");
            dataspace = dataset.getSpace();
            dataspace.selectHyperslab(H5S_SELECT_SET, count,
                                      offset, stride, block);
            dataset.write(temp, NATIVE_REAL_T,
                          memspace, dataspace);

            delete[] temp;
        }

    }

}


// void BunchMonitor::open(){}
void BunchMonitor::init_buffer()
{
    b_np_alive.resize(fBufferTime);
    b_mean_dt.resize(fBufferTime);
    b_mean_dE.resize(fBufferTime);
    b_sigma_dt.resize(fBufferTime);
    b_sigma_dE.resize(fBufferTime);
    b_epsn_rms.resize(fBufferTime);

    if (fGaussian)
        b_bl_gauss.resize(fBufferTime);

    if (fPL != NULL) {
        b_PL_omegaRF.resize(fBufferTime);
        b_PL_phiRF.resize(fBufferTime);
        b_PL_bunch_phase.resize(fBufferTime);
        b_PL_phase_corr.resize(fBufferTime);
        b_PL_omegaRF_corr.resize(fBufferTime);
        b_SL_dphiRF.resize(fBufferTime);
        b_RL_drho.resize(fBufferTime);
    }
    if (fNoiseFB != NULL) {
        b_LHCnoiseFB_factor.resize(fBufferTime);
        b_LHCnoiseFB_bl.resize(fBufferTime);
        if (!fNoiseFB->fBlMeasBBB.empty())
            b_LHCnoiseFB_bl_bbb.resize(fBufferTime);
    }

}


void BunchMonitor::write_buffer()
{

    int i = fITurn % fBufferTime;

    b_np_alive[i] = fBeam->n_macroparticles_alive();
    b_mean_dt[i] = fBeam->mean_dt;
    b_mean_dE[i] = fBeam->mean_dE;
    b_sigma_dt[i] = fBeam->sigma_dt;
    b_sigma_dE[i] = fBeam->sigma_dE;
    b_epsn_rms[i] = fBeam->epsn_rms_l;

    if (fGaussian)
        b_bl_gauss[i] = fSlices->bl_gauss;

    if (fPL != NULL) {
        b_PL_omegaRF[i] = fRfP->omega_RF[0][fITurn];
        b_PL_phiRF[i] = fRfP->phi_RF[0][fITurn];
        b_PL_bunch_phase[i] = fPL->phi_beam;
        b_PL_phase_corr[i] = fPL->dphi;
        b_PL_omegaRF_corr[i] = fPL->domega_RF;
        b_SL_dphiRF[i] = fRfP->dphi_RF[0];
        b_RL_drho[i] = fPL->drho;
    }

    if (fNoiseFB != NULL) {
        b_LHCnoiseFB_factor[i] = fNoiseFB->fX;
        b_LHCnoiseFB_bl[i] = fNoiseFB->fBlMeas;
        if (!fNoiseFB->fBlMeasBBB.empty())
            b_LHCnoiseFB_bl_bbb[i] =
                std::vector<real_t> (fNoiseFB->fBlMeasBBB.begin(),
                                     fNoiseFB->fBlMeasBBB.end());
        // b_LHCnoiseFB_bl_bbb[i] = fNoiseFB->fBlMeasBBB;
    }
}

void BunchMonitor::close()
{
    fH5Group->close();
    fH5File->close();
}
