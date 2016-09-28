#include <blond/plots/plot_llrf.h>
#include <blond/python.h>
#include <blond/monitors/Monitors.h>
#include <memory>


void plot_noise_spectrum(f_vector_t &frequency, f_vector_t &spectrum,
                         int sampling, std::string dirname, int figno)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_noise_spectrum");

    auto pFigNo = python::convert_int(figno);
    auto pSampling = python::convert_int(sampling);
    auto pDirname = python::convert_string(dirname);

    auto pFreq = python::convert_double_array(frequency.data(),
                 frequency.size());

    auto pSpectrum = python::convert_double_array(spectrum.data(),
                     spectrum.size());

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pFreq,
                                            pSpectrum, pSampling,
                                            pDirname, pFigNo, NULL);
    assert(ret);
}


void plot_phase_noise(f_vector_t &time, f_vector_t &dphi, int sampling,
                      std::string dirname, int figno)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_phase_noise");

    auto pFigNo = python::convert_int(figno);
    auto pSampling = python::convert_int(sampling);
    auto pDirname = python::convert_string(dirname);

    auto pTime = python::convert_double_array(time.data(),
                 time.size());

    auto pDPhi = python::convert_double_array(dphi.data(),
                 dphi.size());

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pTime, pDPhi, pSampling,
                                            pDirname, pFigNo, NULL);
    assert(ret);
}


void plot_PL_bunch_phase(RfParameters *RfP, std::string h5data,
                         int output_freq, std::string dirname)
{

    python::import();

    auto pFunc = python::import("plot_llrf", "plot_PL_bunch_phase");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *)read_1D(h5data,
                                  "Beam/PL_bunch_phase",
                                  "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_PL_RF_phase(RfParameters *RfP, std::string h5data,
                      int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_PL_RF_phase");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data, "Beam/PL_phiRF",
                                  "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_PL_phase_corr(RfParameters *RfP, std::string h5data,
                        int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_PL_phase_corr");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data,
                                  "Beam/PL_phase_corr", "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_PL_RF_freq(RfParameters *RfP, std::string h5data,
                     int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_PL_RF_freq");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data, "Beam/PL_omegaRF",
                                  "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_PL_freq_corr(RfParameters *RfP, std::string h5data,
                       int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_PL_freq_corr");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data,
                                  "Beam/PL_omegaRF_corr",
                                  "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_RF_phase_error(RfParameters *RfP, std::string h5data,
                         int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_RF_phase_error");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data, "Beam/SL_dphiRF",
                                  "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_RL_radial_error(RfParameters *RfP, std::string h5data,
                          int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_RL_radial_error");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data, "Beam/RL_drho",
                                  "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_COM_motion(RfParameters *RfP, std::string h5data,
                     int output_freq, std::string dirname)
{

    python::import();

    auto pFunc = python::import("plot_llrf", "plot_COM_motion");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data, "Beam/mean_dt",
                                  "double", &dims));

    auto pMeanDt = python::convert_double_array(array.get(), dims);

    array = std::unique_ptr<double> ((double *) read_1D(h5data,
                                     "Beam/mean_dE", "double", &dims));

    auto pMeanDE = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pMeanDt,
                                            pMeanDE, pOutputFreq, pDirname,
                                            NULL);
    assert(ret);
}


void plot_LHCNoiseFB(RfParameters *RfP, std::string h5data,
                     int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_LHCNoiseFB");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data,
                                  "Beam/LHC_noise_FB_factor",
                                  "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);
}



void plot_LHCNoiseFB_FWHM(RfParameters *RfP, std::string h5data,
                          int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_LHCNoiseFB_FWHM");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> array((double *) read_1D(h5data,
                                  "Beam/LHC_noise_FB_bl", "double", &dims));

    auto ph5Data = python::convert_double_array(array.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);
}


// TODO fix this to read & convert 2d array
void plot_LHCNoiseFB_FWHM_bbb(RfParameters *RfP, std::string h5data,
                              int output_freq, std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_llrf", "plot_LHCNoiseFB_FWHM_bbb");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims[2];
    std::unique_ptr<double> array((double *) read_2D(h5data,
                                  "Beam/LHC_noise_FB_bl_bbb",
                                  "double", dims));

    f_vector_2d_t temp;
    for (uint i = 0; i < dims[0]; i++) {
        temp.push_back(f_vector_t(&array.get()[i * dims[1]],
                                  &array.get()[(i + 1) * dims[1]]));
    }

    auto ph5Data = python::convert_double_2d_array(temp);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, ph5Data,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);
}
