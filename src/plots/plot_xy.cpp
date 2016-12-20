/*
 * plot_xy.cpp
 *
 *  Created on: Dec 20, 2016
 *      Author: alasheen
 */


#include <blond/plots/plot_xy.h>
#include <blond/python.h>
#include <blond/utilities.h>

int plot_xy(f_vector_t xArray, f_vector_t yArray, std::string figureName, std::string dirName, std::string saveFigName, double xMax) {

  
  python::import();
  auto pFunc = python::import("plot_xy", "plot_xy");

  auto py_xArray = python::convert_double_array(xArray.data(), xArray.size());
  auto py_yArray = python::convert_double_array(yArray.data(), yArray.size());
  auto py_figureName = python::convert_string(figureName);
  auto py_dirName = python::convert_string(dirName);
  auto py_saveFigName = python::convert_string(saveFigName);
  auto py_xMax = python::convert_double(xMax);

  auto ret = PyObject_CallFunctionObjArgs(pFunc, py_xArray, py_yArray, py_figureName, py_dirName, py_saveFigName, py_xMax, NULL);
  return ret != nullptr;

}

