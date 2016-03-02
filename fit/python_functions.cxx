////////////////////////////////////////
// This file contains python_functions implementation
////////////////////////////////////////
#include <Python.h>
#include "fit_factory.h"

using namespace std;
static fit_factory* ff;
static fit_exponential* expo;
static fit_polynomial* poly;
static fit_gaussian * gaussian;
static fit_radial_basis_function* rbf;
static fit_lm_solver *lm;
static fit_amoeba_solver *amoeba;
static fit_svd_solver *svd;

PyObject* fit(PyObject* /*self*/, PyObject *args) {
	PyObject* xList_obj, *yList_obj;
	PyObject* listItemX_obj, *listItemY_obj;
  char* fnname, *fsname;
  if (!PyArg_ParseTuple(args, "O!O!ss",
    &PyList_Type, &xList_obj, &PyList_Type, &yList_obj, &fnname, &fsname))
    return NULL;

	unsigned num_elem = PyList_Size(xList_obj);	
  if (num_elem <= 0) return NULL;
	vnl_vector<double> xi(num_elem), y(num_elem);
	for (unsigned i = 0; i < num_elem; i++) { 	/* iterate over items of the list, grabbing PyObjects, and converting to doubles */
		listItemX_obj = PyList_GetItem(xList_obj, i);
		xi[i] = (PyFloat_AsDouble(listItemX_obj)); /* convert to double*/
		listItemY_obj = PyList_GetItem(yList_obj, i);
		y[i] = PyFloat_AsDouble(listItemY_obj);
	}

  vnl_vector<double> fit_params;
  ff->fit(xi, y, fnname, fsname, fit_params);
  PyObject* params = PyList_New(fit_params.size());  //parameters
  for (unsigned i = 0; i < fit_params.size(); i++) {
    PyList_SetItem(params, i, Py_BuildValue("d", fit_params[i]));
  }

	PyObject* y_vals = PyList_New(num_elem); // fitted model y_vals
  for (unsigned i = 0; i < num_elem; i++) {
		PyList_SetItem(y_vals, i, Py_BuildValue("d", ff->fn()->g(xi[i])));
	}

	PyObject* outlist = PyTuple_New(2); // return a tuple: (parameters, y_vals)
	PyTuple_SetItem(outlist, 0, Py_BuildValue("O", params));
	PyTuple_SetItem(outlist, 1, Py_BuildValue("O", y_vals));
	return outlist;
}


PyObject* register_fit_factory(PyObject* /*self*/) {
  ff = new fit_factory;	
  expo = new fit_exponential;	ff->register_function(expo);
  poly = new fit_polynomial;	ff->register_function(poly);
  gaussian = new fit_gaussian; ff->register_function(gaussian);
  rbf = new fit_radial_basis_function; ff->register_function(rbf);
  lm = new fit_lm_solver;	ff->register_solver(lm);
  amoeba = new fit_amoeba_solver;	ff->register_solver(amoeba);
  svd = new fit_svd_solver; ff->register_solver(svd);
  Py_RETURN_NONE;
}

PyObject* remove_fit_factory(PyObject* /*self*/) {
  delete lm; delete amoeba; delete svd;
  delete gaussian; delete expo; delete poly; delete rbf;
  delete ff;
  Py_RETURN_NONE;
}

static PyMethodDef py_fit_methods[] = {		
	{ "register_fit_factory", (PyCFunction)register_fit_factory, METH_NOARGS, "Register Fit Factory" },
	{ "remove_fit_factory", (PyCFunction)remove_fit_factory, METH_NOARGS, "Remove Fit Factory" },
  { "fit", (PyCFunction)fit, METH_VARARGS, "fit to gaussian function" },
	{ NULL, NULL, 0, NULL },
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef py_fit = { PyModuleDef_HEAD_INIT, "py_fit", NULL, -1, py_fit_methods};
PyMODINIT_FUNC PyInit_py_fit(void) {
  return PyModule_Create(&py_fit);
}

#else // Python 2
PyMODINIT_FUNC initpy_fit(void) {
  (void)Py_InitModule("py_fit", py_fit_methods);
}
#endif
