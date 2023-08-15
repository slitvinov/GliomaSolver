#define PY_SSIZE_T_CLEAN
#include "lib.h"
#include <Python.h>

#define max(a, b) ((a) > (b) ? (a) : (b))

static PyObject *run(PyObject *self, PyObject *args) {
  Py_buffer gm_view, wm_view, hg_view;
  Py_ssize_t indices[3] = {0, 0, 0};
  double whenToWrite, whenToWriteOffset, tend, h, t, L;
  int step, gpd;
  PyObject *GM, *WM, *HG;
  float *hg;
  struct Brain *brain;
  struct BrainParams params;

  if (!PyArg_ParseTuple(args, "OOO", &GM, &WM, &HG))
    return NULL;
  if (!PyObject_CheckBuffer(GM) || !PyObject_CheckBuffer(WM) ||
      !PyObject_CheckBuffer(HG))
    return NULL;
  if (PyObject_GetBuffer(GM, &gm_view, PyBUF_F_CONTIGUOUS) == -1)
    return NULL;
  if (PyObject_GetBuffer(WM, &wm_view, PyBUF_F_CONTIGUOUS) == -1)
    return NULL;
  if (PyObject_GetBuffer(HG, &hg_view, PyBUF_F_CONTIGUOUS | PyBUF_WRITABLE) ==
      -1)
    return NULL;
  if (gm_view.ndim != 3 || wm_view.ndim != 3 || hg_view.ndim != 3) {
    PyErr_SetString(PyExc_ValueError, "GM, WM, HG: ndim != 3");
    return NULL;
  }
  if (gm_view.shape[0] != wm_view.shape[0] ||
      gm_view.shape[1] != wm_view.shape[1] ||
      gm_view.shape[2] != wm_view.shape[2]) {
    PyErr_SetString(PyExc_ValueError, "GM, WM: dimensions do not match");
    return NULL;
  }

  if (gm_view.itemsize != sizeof(float) || wm_view.itemsize != sizeof(float) ||
      hg_view.itemsize != sizeof(float)) {
    PyErr_SetString(PyExc_ValueError, "GM, WM, HG: wrong type");
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

  params.blocksPerDimension = 16;
  params.ic[0] = 0.6497946102507519;
  params.ic[1] = 0.5908331665234543;
  params.ic[2] = 0.3715947899171972;
  gpd = (_BLOCKSIZE_)*params.blocksPerDimension;

  if (hg_view.shape[0] != gpd || hg_view.shape[1] != gpd ||
      hg_view.shape[2] != gpd) {
    PyErr_Format(PyExc_ValueError,
                 "HG: dimensions do not match, expecting %dx%dx%d", gpd, gpd,
                 gpd);
    return NULL;
  }

  h = 1. / gpd;
  L = max(params.n[0], max(params.n[1], params.n[2])) * 0.1;
  params.Dw = 0.0013 / (L * L);
  params.rho = 0.025;
  tend = 300;
  params.dt = 0.99 * h * h / (2. * 3 * params.Dw);
  brain_ini(&params, &brain);
  t = 0.0;
  step = 0;
  while (t <= tend) {
    if (PyErr_CheckSignals() != 0)
      break;
    brain_step(brain);
    t += params.dt;
    step++;
  }
  hg = PyBuffer_GetPointer(&hg_view, indices);
  brain_project(brain, hg);
  brain_fin(brain);
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
