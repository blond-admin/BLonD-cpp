#include <blond/plots/plot_beams.h>
#include <blond/python.h>


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

    auto pRfPCounter = python::convert_int(RfP->counter);
    auto pRfPOmegaRf = python::convert_double(RfP->omega_RF[0][RfP->counter]);
    auto pRfPPhiRf = python::convert_double(RfP->phi_RF[0][RfP->counter]);
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

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pRfPCounter, pRfPOmegaRf,
                                            pRfPPhiRf, pBeamId, pBeamDt,
                                            pBeamDE, pXMin, pXMax, pYMin, pYMax,
                                            pXUnit, pSampling, pSeparatrixPlot,
                                            pHistogramsPlot, pDirname, pAlpha,
                                            NULL);
    assert(ret);

    // python::finalize();

}