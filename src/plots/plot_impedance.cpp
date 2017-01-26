#include <blond/plots/plot_impedance.h>
#include <blond/python.h>
using namespace blond;

int blond::plot_impedance_vs_frequency(int counter,
                                InducedVoltageFreq *indVoltFreq,
                                Slices *slices,
                                std::string option1,
                                std::string option2,
                                std::string option3,
                                std::string style,
                                f_vector_t cut_left_right,
                                f_vector_t cut_up_down,
                                std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_impedance",
                                "plot_impedance_vs_frequency");

    auto pCounter = python::convert_int(counter);
    auto pDirname = python::convert_string(dirname);
    auto pStyle = python::convert_string(style);
    auto pOption1 = python::convert_string(option1);
    auto pOption2 = python::convert_string(option2);
    auto pOption3 = python::convert_string(option3);
    auto pCutLR = python::convert_double_array(cut_left_right.data(),
                  cut_left_right.size());

    auto pCutUD = python::convert_double_array(cut_up_down.data(),
                  cut_up_down.size());

    auto pTotalImpedance = python::convert_complex_array(
                               indVoltFreq->fTotalImpedance.data(),
                               indVoltFreq->fTotalImpedance.size());

    auto pBeamSpec = python::convert_complex_array(
                         slices->fBeamSpectrum.data(),
                         slices->fBeamSpectrum.size());

    auto pBeamSpecFreq = python::convert_double_array(
                             slices->fBeamSpectrumFreq.data(),
                             slices->fBeamSpectrumFreq.size());

    auto pFreqArray = python::convert_double_array(
                          indVoltFreq->fFreqArray.data(),
                          indVoltFreq->fFreqArray.size());

    f_vector_2d_t freq_array_2d;
    f_vector_2d_t impedance_real_2d;
    f_vector_2d_t impedance_imag_2d;

    for (auto &intensity : indVoltFreq->fImpedanceSourceList) {
        auto p = dynamic_cast<InputTable *>(intensity);
        if (p != nullptr && option3 == "freq_table") {
            freq_array_2d.push_back(p->fFrequencyArrayLoaded);
            impedance_real_2d.push_back(p->fReZArrayLoaded);
            impedance_imag_2d.push_back(p->fImZArrayLoaded);
        } else {

            freq_array_2d.push_back(indVoltFreq->fFreqArray);
            f_vector_t real, imag;
            for (const auto &z : intensity->fImpedance) {
                real.push_back(z.real());
                imag.push_back(z.imag());
            }
            impedance_real_2d.push_back(real);
            impedance_imag_2d.push_back(imag);
        }
    }

    auto pFreqArray2D = python::convert_double_2d_array(freq_array_2d);
    auto pIMpedanceReal2D = python::convert_double_2d_array(impedance_real_2d);
    auto pIMpedanceImag2D = python::convert_double_2d_array(impedance_imag_2d);


    auto ret = PyObject_CallFunctionObjArgs(pFunc, pCounter, pFreqArray,
                                            pTotalImpedance, pBeamSpec,
                                            pBeamSpecFreq, pFreqArray2D,
                                            pIMpedanceReal2D, pIMpedanceImag2D,
                                            pOption1, pOption2, pOption3,
                                            pStyle, pCutLR, pCutUD,
                                            pDirname, NULL);
    return ret != nullptr;

}



int blond::plot_induced_voltage_vs_bin_centers(int counter,
                                        TotalInducedVoltage *totIndVolt,
                                        Slices *slices,
                                        std::string style ,
                                        std::string dirname)
{
    python::import();

    auto pFunc = python::import("plot_impedance",
                                "plot_induced_voltage_vs_bin_centers");

    auto pCounter = python::convert_int(counter);
    auto pDirname = python::convert_string(dirname);
    auto pStyle = python::convert_string(style);

    auto pBinCenters = python::convert_double_array(
                           slices->bin_centers.data(),
                           slices->bin_centers.size());

    auto pInducedVoltage = python::convert_double_array(
                               totIndVolt->fInducedVoltage.data(),
                               totIndVolt->fInducedVoltage.size());

    auto ret = PyObject_CallFunctionObjArgs(pFunc, pCounter, pBinCenters,
                                            pInducedVoltage, pStyle,
                                            pDirname, NULL);
    return ret != nullptr;
}
