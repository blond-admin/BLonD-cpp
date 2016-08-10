
#ifndef INCLUDE_BLOND_PYTHON_H_
#define INCLUDE_BLOND_PYTHON_H_

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/core/include/numpy/arrayobject.h>
#include <iostream>
#include <map>


namespace python {

    // static bool 

    static inline int initialize()
    {
        Py_Initialize();
        // import_array1(0);
        return 0;
    }

    static inline int import()
    {
        // Py_Initialize();
        import_array1(0);
        return 0;
    }

    static inline void finalize()
    {
        Py_Finalize();
    }

    static inline PyObject *import(std::string module, std::string function)
    {
        auto pModule = PyImport_ImportModule(module.c_str());
        assert(pModule);

        auto pFunc = PyObject_GetAttrString(pModule, function.c_str());
        assert(pFunc);

        return pFunc;
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
        auto pVar = (PyArrayObject *) PyArray_FromDimsAndData(1, dims, NPY_INT,
                    reinterpret_cast<char *>(array));

        assert(pVar);
        return pVar;
    }

    static inline PyArrayObject *convert_double_array(double *array, int size)
    {
        int dims[1] = {size};
        auto pVar = (PyArrayObject *) PyArray_FromDimsAndData(1, dims, NPY_DOUBLE,
                    reinterpret_cast<char *>(array));

        assert(pVar);
        return pVar;
    }

    static inline PyObject *convert_dictionary(std::map<std::string, std::string> map)
    {
        // int dims[1] = {size};
        // auto pArray = (PyArrayObject *) PyArray_FromDimsAndData(1, dims, NPY_DOUBLE,
        //               reinterpret_cast<char *>(array));
        if (map.size() > 0) {
            auto pDict = PyDict_New();
            assert(pDict);
            for (const auto &pair : map) {
                auto pKey = PyString_FromString(pair.first.c_str());
                auto pVal = PyString_FromString(pair.second.c_str());
                assert(PyDict_SetItem(pDict, pKey, pVal) == 0);
            }
            return pDict;
        } else {
            return Py_None;
        }
    }

}


#endif /* INCLUDE_BLOND_PYTHON_H_ */
