#include <blond/plots/plot_beams.h>
#include <blond/python.h>
#include <blond/monitors/Monitors.h>

void plot_long_phase_space(GeneralParameters *GP, RfParameters *RfP,
                           Beams *Beam, double xmin, double xmax,
                           double ymin, double ymax, std::string xunit,
                           int sampling, bool separatrix_plot,
                           bool histograms_plot, std::string dirname,
                           int alpha)
{
    // python::initialize();
    python::import();

    auto pFunc = python::import("plot_beams", "plot_long_phase_space");

    uint turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pRfPOmegaRf0 = python::convert_double(RfP->omega_RF[0][turn]);
    auto pRfPPhiRf0 = python::convert_double(RfP->phi_RF[0][turn]);
    auto pBeamDE = python::convert_double_array(Beam->dE.data(), Beam->dE.size());
    auto pBeamDt = python::convert_double_array(Beam->dt.data(), Beam->dt.size());
    auto pBeamId = python::convert_int_array(Beam->id.data(), Beam->id.size());

    auto pXMin = python::convert_double(xmin);
    auto pXMax = python::convert_double(xmax);
    auto pYMin = python::convert_double(ymin);
    auto pYMax = python::convert_double(ymax);
    auto pXUnit = python::convert_string(xunit);
    auto pSampling = python::convert_int(sampling);
    auto pSeparatrixPlot = python::convert_bool(separatrix_plot);
    auto pHistogramsPlot = python::convert_bool(histograms_plot);
    auto pDirname = python::convert_string(dirname);
    auto pAlpha = python::convert_int(alpha);

    if (separatrix_plot == true) {
        auto pGPNSections = python::convert_int(GP->n_sections);
        auto pGPCharge = python::convert_double(GP->charge);
        auto pGPTRev = python::convert_double(GP->t_rev[turn]);

        f_vector_t temp;
        for (auto &row : RfP->voltage) temp.push_back(row[turn]);
        auto pRfPVoltage = python::convert_double_array(temp.data(), temp.size());

        temp.clear();
        for (auto &row : RfP->omega_RF) temp.push_back(row[turn]);
        auto pRfPOmegaRf = python::convert_double_array(temp.data(), temp.size());

        temp.clear();
        for (auto &row : RfP->phi_RF) temp.push_back(row[turn]);
        auto pRfPPhiRf = python::convert_double_array(temp.data(), temp.size());

        auto pRfPEta0 = python::convert_double(RfP->eta_0(turn));
        auto pRfPBeta = python::convert_double(RfP->beta(turn));
        auto pRfPEnergy = python::convert_double(RfP->energy(turn));
        auto pRfPNRf = python::convert_int(RfP->n_rf);
        auto pRfPHarmonic = python::convert_double(RfP->harmonic[0][turn]);
        auto pRfPPhiS = python::convert_double(RfP->phi_s[turn]);
        auto pRfPEIncrement = turn >= RfP->E_increment.size() ?
                              python::convert_double(RfP->E_increment.back()) :
                              python::convert_double(RfP->E_increment[turn]);
        auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pRfPOmegaRf0,
                                                pRfPPhiRf0, pBeamId, pBeamDt,
                                                pBeamDE, pXMin, pXMax, pYMin, pYMax,
                                                pXUnit, pSampling, pSeparatrixPlot,
                                                pGPNSections, pGPCharge, pGPTRev,
                                                pRfPVoltage, pRfPOmegaRf, pRfPPhiRf,
                                                pRfPEta0, pRfPBeta, pRfPEnergy, pRfPNRf,
                                                pRfPHarmonic, pRfPPhiS, pRfPEIncrement,
                                                pHistogramsPlot, pDirname, pAlpha,
                                                NULL);
        assert(ret);

    } else {
        auto pGPNSections = python::get_none();
        auto pGPCharge = python::get_none();
        auto pGPTRev = python::get_none();
        auto pRfPVoltage = python::get_none();
        auto pRfPOmegaRf = python::get_none();
        auto pRfPPhiRf = python::get_none();
        auto pRfPEta0 = python::get_none();
        auto pRfPBeta = python::get_none();
        auto pRfPEnergy = python::get_none();
        auto pRfPNRf = python::get_none();
        auto pRfPHarmonic = python::get_none();
        auto pRfPPhiS = python::get_none();
        auto pRfPEIncrement = python::get_none();
        auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pRfPOmegaRf0,
                                                pRfPPhiRf0, pBeamId, pBeamDt,
                                                pBeamDE, pXMin, pXMax, pYMin, pYMax,
                                                pXUnit, pSampling, pSeparatrixPlot,
                                                pGPNSections, pGPCharge, pGPTRev,
                                                pRfPVoltage, pRfPOmegaRf, pRfPPhiRf,
                                                pRfPEta0, pRfPBeta, pRfPEnergy, pRfPNRf,
                                                pRfPHarmonic, pRfPPhiS, pRfPEIncrement,
                                                pHistogramsPlot, pDirname, pAlpha,
                                                NULL);
        assert(ret);
    }



    // python::finalize();

}


void plot_bunch_length_evol(RfParameters *RfP, std::string h5data,
                            int output_freq, std::string dirname)
{

    python::import();

    auto pFunc = python::import("plot_beams", "plot_bunch_length_evol");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> sigma_dt((double *) read_1D(h5data, "Beam/sigma_dt",
                                     "double", &dims));

    auto pSigmaDt = python::convert_double_array(sigma_dt.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pSigmaDt,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}

// NOTE removed unused variable Slice
void plot_bunch_length_evol_gaussian(RfParameters *RfP,
                                     std::string h5data, int output_freq,
                                     std::string dirname)
{

    python::import();

    auto pFunc = python::import("plot_beams", "plot_bunch_length_evol_gaussian");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);

    hsize_t dims;
    std::unique_ptr<double> bl_gauss((double *) read_1D(h5data,
                                     "Beam/bunch_length_gaussian",
                                     "double", &dims));

    auto pBlGauss = python::convert_double_array(bl_gauss.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pBlGauss,
                                            pOutputFreq, pDirname, NULL);
    assert(ret);

}


void plot_position_evol(RfParameters *RfP, std::string h5data,
                        int output_freq, std::string style,
                        std::string dirname)
{

    python::import();

    auto pFunc = python::import("plot_beams", "plot_position_evol");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);
    auto pStyle = python::convert_string(style);

    hsize_t dims;
    std::unique_ptr<double> mean_dt((double *) read_1D(h5data, "Beam/mean_dt",
                                    "double", &dims));

    auto pMeanDt = python::convert_double_array(mean_dt.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pMeanDt,
                                            pOutputFreq, pStyle, pDirname,
                                            NULL);
    assert(ret);
}


void plot_energy_evol(RfParameters *RfP, std::string h5data,
                      int output_freq, std::string style,
                      std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_beams", "plot_energy_evol");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);
    auto pStyle = python::convert_string(style);

    hsize_t dims;
    std::unique_ptr<double> mean_dE((double *) read_1D(h5data, "Beam/mean_dE",
                                    "double", &dims));

    auto pMeanDE = python::convert_double_array(mean_dE.get(), dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pMeanDE,
                                            pOutputFreq, pStyle, pDirname,
                                            NULL);
    assert(ret);
}

void plot_transmitted_particles(RfParameters *RfP, std::string h5data,
                                int output_freq, std::string style,
                                std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_beams", "plot_transmitted_particles");

    int turn = RfP->counter;
    auto pRfPCounter = python::convert_int(turn);
    auto pOutputFreq = python::convert_int(output_freq);
    auto pDirname = python::convert_string(dirname);
    auto pStyle = python::convert_string(style);

    hsize_t dims;
    std::unique_ptr<int> n_macroparticles_alive((int *) read_1D(h5data,
            "Beam/n_macroparticles_alive",
            "int", &dims));

    auto pPartsAlive = python::convert_int_array(n_macroparticles_alive.get(),
                       dims);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pPartsAlive,
                                            pOutputFreq, pStyle, pDirname,
                                            NULL);
    assert(ret);
}
