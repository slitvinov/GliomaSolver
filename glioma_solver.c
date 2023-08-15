#define PY_SSIZE_T_CLEAN
#include "lib.h"
#include <Python.h>

static PyObject *run(PyObject *self, PyObject *args) {
  struct BrainParams params;
  PyObject *GM, *WM;
  const Py_buffer gm_view, wm_view;
  const Py_ssize_t indices[3] = {0, 0, 0};

  if (!PyArg_ParseTuple(args, "OO", &GM, &WM))
    return NULL;
  if (!PyObject_CheckBuffer(GM))
    return NULL;
  if (!PyObject_CheckBuffer(WM))
    return NULL;
  if (PyObject_GetBuffer(GM, &gm_view, PyBUF_F_CONTIGUOUS) == -1)
    return NULL;
  if (PyObject_GetBuffer(WM, &wm_view, PyBUF_F_CONTIGUOUS) == -1)
    return NULL;
  if (gm_view.ndim != 3) {
    PyErr_SetString(PyExc_ValueError, "GM: ndim != 3");
    return NULL;
  }
  if (wm_view.ndim != 3) {
    PyErr_SetString(PyExc_ValueError, "WM: ndim != 3");
    return NULL;
  }
  if (gm_view.shape[0] != wm_view.shape[0] ||
      gm_view.shape[1] != wm_view.shape[1] ||
      gm_view.shape[2] != wm_view.shape[2]) {
    PyErr_SetString(PyExc_ValueError, "WM an GM: dimensions do not match");
    return NULL;
  }

  if (gm_view.itemsize != sizeof(float) || wm_view.itemsize != sizeof(float)) {
    PyErr_SetString(PyExc_ValueError, "WM an GM: wrong type");
    return NULL;
  }

  params.GM = PyBuffer_GetPointer(&gm_view, indices);
  params.WM = PyBuffer_GetPointer(&gm_view, indices);
  params.n[0] = gm_view.shape[0];
  params.n[1] = gm_view.shape[1];
  params.n[2] = gm_view.shape[2];

  fprintf(stderr, "%d %d %d\n", params.n[0], params.n[1], params.n[2]);
  fprintf(stderr, "params.GM: %g %g\n", params.GM[0],
          params.GM[256 * 256 * 256 - 1]);

  return PyFloat_FromDouble(0.0);
}

static PyMethodDef Glioma_SolverMethods[] = {
    {"run", run, METH_VARARGS, PyDoc_STR("Run simulations")},
    {NULL, NULL, 0, NULL},
};

static PyModuleDef glioma_solvermodule = {
    PyModuleDef_HEAD_INIT, "glioma_solver", PyDoc_STR("Glioma Solver"), -1,
    Glioma_SolverMethods,
};

PyMODINIT_FUNC PyInit_glioma_solver(void) {
  return PyModule_Create(&glioma_solvermodule);
}
