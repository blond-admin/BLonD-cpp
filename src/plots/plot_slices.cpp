#include <blond/plots/plot_slices.h>
#include <blond/python.h>


void plot_beam_profile(Slices *Slices,
                       int counter,
                       std::string style,
                       std::string dirname)
{
    // python::initialize();
    python::import();
    auto pFunc = python::import("plot_slices", "plot_beam_profile");

    auto pBinCenters = python::convert_double_array(Slices->bin_centers.data(),
                       Slices->bin_centers.size());
    auto pNMacropaticles = python::convert_int_array(Slices->n_macroparticles.data(),
                           Slices->n_macroparticles.size());
    auto pCounter = python::convert_int(counter);
    auto pStyle = python::convert_string(style);
    auto pDirname = python::convert_string(dirname);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pBinCenters, pNMacropaticles,
                                            pCounter, pStyle,  pDirname,
                                            NULL);
    assert(ret);
    // python::finalize();
}


// TODO: Implement beam_profile_derivative
void plot_beam_profile_derivative(Slices *Slices,
                                  int counter,
                                  std::string style,
                                  std::string dirname,
                                  int_vector_t numbers)
{
    // python::import();
    // auto pFunc = python::import("plot_slices", "plot_beam_profile_derivative");

    // auto pBinCenters = python::convert_double_array(Slices->bin_centers.data(),
    //                    Slices->bin_centers.size());
    // auto pNMacropaticles = python::convert_int_array(Slices->n_macroparticles.data(),
    //                        Slices->n_macroparticles.size());
    // auto pNumbers = python::convert_int_array(numbers.data(), numbers.size());
    // auto pCounter = python::convert_int(counter);
    // auto pStyle = python::convert_string(style);
    // auto pDirname = python::convert_string(dirname);

    // auto ret = PyObject_CallFunctionObjArgs(pFunc, pBinCenters, pNMacropaticles,
    //                                         pDerivative, pX, pCounter, pStyle,
    //                                         pDirname, pNumbers, NULL);
    // assert(ret);
}

void plot_beam_spectrum(Slices *Slices, int counter,
                        std::string style,
                        std::string dirname)
{

    python::import();
    auto pFunc = python::import("plot_slices", "plot_beam_spectrum");

    auto pBeamSpectrumFreq = python::convert_double_array(
                                 Slices->fBeamSpectrumFreq.data(),
                                 Slices->fBeamSpectrumFreq.size());
    auto pBeamSpectrum = python::convert_complex_array(
                             Slices->fBeamSpectrum.data(),
                             Slices->fBeamSpectrum.size());


    auto pCounter = python::convert_int(counter);
    auto pStyle = python::convert_string(style);
    auto pDirname = python::convert_string(dirname);

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pBeamSpectrumFreq,
                                            pBeamSpectrum, pCounter, pStyle,
                                            pDirname, NULL);
    assert(ret);

}
