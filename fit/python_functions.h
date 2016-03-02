////////////////////////////////////////
// This file contains python_functions class
////////////////////////////////////////
#ifndef python_functions_h_
#define python_functions_h_
#include <iostream>
#include <Python.h>

PyObject* fit(PyObject* /*self*/, PyObject *args);
PyObject* register_fit_factory(PyObject* /*self*/);
PyObject* remove_fit_factory(PyObject* /*self*/);

#endif //python_functions_h_