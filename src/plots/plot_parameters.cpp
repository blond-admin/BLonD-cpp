
#include <blond/plots/plot_parameters.h>
#include <blond/python.h>


void plot_voltage_programme(f_vector_t &time, f_vector_t &voltage,
                            int sampling, std::string dirname, int figno)
{
    // python::initialize();
    python::import();

    auto pFunc = python::import("plot_parameters", "plot_voltage_programme");

    auto pTime = python::convert_double_array(time.data(), time.size());
    auto pVoltage = python::convert_double_array(voltage.data(), voltage.size());
    auto pSampling = python::convert_int(sampling);
    auto pDirname = python::convert_string(dirname);
    auto pFigno = python::convert_int(figno);


    auto ret = PyObject_CallFunctionObjArgs(pFunc, pTime, pVoltage, pSampling,
                                            pDirname, pFigno, NULL);
    assert(ret);



}