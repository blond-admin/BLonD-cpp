
#ifndef INCLUDE_BLOND_PYTHON_H_
#define INCLUDE_BLOND_PYTHON_H_

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#ifdef __GNUC__
// Avoid tons of warnings with root code
#pragma GCC system_header
#endif

#include <Python.h>
#include <numpy/core/include/numpy/arrayobject.h>
#include <map>
#include <blond/configuration.h>
#include <cstdlib>

namespace python {

    static inline int initialize()
    {
        putenv("PYTHONPATH=" PYPATH);
        Py_SetPythonHome(PYHOME);
        Py_Initialize();
        return 0;
    }

    static inline int import()
    {
        import_array1(0);
        return 0;
    }

    static inline void finalize()
    {
        Py_Finalize();
    }

    static inline PyObject *import(std::string module, std::string function)
    {
        // TODO try to prepend the path to python folder
        auto pModule = PyImport_ImportModule(module.c_str());
        assert(pModule);

        auto pFunc = PyObject_GetAttrString(pModule, function.c_str());
        assert(pFunc);

        return pFunc;
    }

    static inline PyObject *get_none()
    {
        return Py_None;
    }

    static inline PyObject *convert_double(double value)
    {
        auto pVar = PyFloat_FromDouble(value);
        assert(pVar);
        return pVar;
    }


    static inline PyObject *convert_int(int value)
    {
        auto pVar = PyInt_FromLong(value);
        assert(pVar);
        return pVar;
    }

    static inline PyObject *convert_string(std::string &value)
    {
        auto pVar = PyString_FromString(value.c_str());
        assert(pVar);
        return pVar;
    }

    static inline PyObject *convert_bool(bool value)
    {
        auto pVar = PyBool_FromLong((long) value);
        assert(pVar);
        return pVar;
    }


    static inline PyArrayObject *convert_int_array(int *array, int size)
    {
        int dims[1] = {size};
        auto pVar = (PyArrayObject *) PyArray_FromDimsAndData(1, dims,
                    NPY_INT, (char *)array);

        assert(pVar);
        return pVar;
    }

    static inline PyArrayObject *convert_double_array(double *array, int size)
    {
        int dims[1] = {size};
        auto pVar = (PyArrayObject *) PyArray_FromDimsAndData(1, dims,
                    NPY_DOUBLE, (char *)array);

        assert(pVar);
        return pVar;
    }

    static inline PyObject *convert_string_array(std::string *array, int size)
    {
        std::string result = array[0];
        for (int i = 1; i < size; i++)
            result += " " + array[i];
        auto pVar = PyString_FromString(result.c_str());
        assert(pVar);
        return pVar;
    }


    /*
    static inline PyArrayObject *convert_double_2d_array(f_vector_2d_t &v)
    {
        int dims[2] = {v.size(), v.front().size()};
        auto array = new double[dims[0] * dims[1]];
        int count = 0;
        for (const auto &row : v) {
            assert(row.size() == dims[1]);
            std::copy(row.begin(), row.end(), &array[count]);
            count += dims[1];
        }
        auto pVar = (PyArrayObject *) PyArray_FromDimsAndData(2, dims, NPY_DOUBLE,
                    (char *)array);
        assert(pVar);
        return pVar;
    }
    */

    static inline PyObject *convert_double_2d_array(f_vector_2d_t &v)
    {
        int xsize = v.size();
        auto pList = PyList_New(xsize);
        assert(pList);
        for (int i = 0; i < xsize; i++) {
            int ysize = v[i].size();
            auto pRow = PyList_New(ysize);
            for (int j = 0; j < ysize; j++) {
                auto pVar = PyFloat_FromDouble(v[i][j]);
                assert(pVar);
                PyList_SET_ITEM(pRow, j, pVar);
            }
            PyList_SET_ITEM(pList, i, pRow);
            // int dims[1] = {v[i].size()};
            // auto pVar = (PyObject *) PyArray_FromDimsAndData(1, dims,
            //             NPY_DOUBLE,
            //             (char *)v[i].data());
            // assert(pVar);
            // PyList_SET_ITEM(pList, i, pVar);
        }
        return pList;
    }



    static inline PyObject *convert_complex_array(std::complex<double> *array,
            int size)
    {
        auto pList = PyList_New(size);
        assert(pList);
        for (int i = 0; i < size; i++) {
            auto z = array[i];
            auto pZ = PyComplex_FromDoubles(z.real(), z.imag());
            assert(pZ);
            PyList_SET_ITEM(pList, i, pZ);
        }
        return pList;
    }

    static inline PyObject *convert_dictionary(std::map<std::string,
            std::string> map)
    {
        if (map.size() > 0) {
            auto pDict = PyDict_New();
            assert(pDict);
            for (const auto &pair : map) {
                auto pKey = PyString_FromString(pair.first.c_str());
                auto pVal = PyString_FromString(pair.second.c_str());
                PyDict_SetItem(pDict, pKey, pVal);
            }
            return pDict;
        } else {
            return Py_None;
        }
    }

}


#endif /* INCLUDE_BLOND_PYTHON_H_ */
