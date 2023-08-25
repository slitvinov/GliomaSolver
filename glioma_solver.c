#define PY_SSIZE_T_CLEAN
#include "lib.h"
#include <Python.h>

#define max(a, b) ((a) > (b) ? (a) : (b))

static PyObject *likelihood(PyObject *self, PyObject *args) {
  Py_ssize_t indices[3] = {0, 0, 0};
  float *model;
  PyObject *MODEL, *PET, *T1c, *FLAIR;
  Py_buffer model_view, pet_view, t1c_view, flair_view;
  struct LikelihoodParams params;
  double ans;

  if (!PyArg_ParseTuple(args, "OOOOddddd", &MODEL, &PET, &T1c, &FLAIR,
                        &params.PETsigma2, &params.PETscale, &params.slope,
                        &params.T1uc, &params.T2uc))
    return NULL;
  if (!PyObject_CheckBuffer(MODEL) || !PyObject_CheckBuffer(PET) ||
      !PyObject_CheckBuffer(T1c) || !PyObject_CheckBuffer(FLAIR))
    return NULL;
  if (PyObject_GetBuffer(MODEL, &model_view, PyBUF_F_CONTIGUOUS) == -1 ||
      PyObject_GetBuffer(PET, &pet_view, PyBUF_F_CONTIGUOUS) == -1 ||
      PyObject_GetBuffer(T1c, &t1c_view, PyBUF_F_CONTIGUOUS) == -1 ||
      PyObject_GetBuffer(FLAIR, &flair_view, PyBUF_F_CONTIGUOUS) == -1)
    return NULL;
  if (model_view.ndim != 3 || pet_view.ndim != 3 || t1c_view.ndim != 3 ||
      flair_view.ndim != 3) {
    PyErr_SetString(PyExc_ValueError, "MODEL, PET, T1c, or FLAIR: ndim != 3");
    return NULL;
  }
  if (model_view.shape[0] != pet_view.shape[0] ||
      model_view.shape[1] != pet_view.shape[1] ||
      model_view.shape[2] != pet_view.shape[2]) {
    PyErr_SetString(PyExc_ValueError, "MODEL, PET: dimensions do not match");
    return NULL;
  }
  if (model_view.shape[0] != t1c_view.shape[0] ||
      model_view.shape[1] != t1c_view.shape[1] ||
      model_view.shape[2] != t1c_view.shape[2]) {
    PyErr_SetString(PyExc_ValueError, "MODEL, T1c: dimensions do not match");
    return NULL;
  }
  if (model_view.shape[0] != flair_view.shape[0] ||
      model_view.shape[1] != flair_view.shape[1] ||
      model_view.shape[2] != flair_view.shape[2]) {
    PyErr_SetString(PyExc_ValueError, "MODEL, FLAIR: dimensions do not match");
    return NULL;
  }
  if (model_view.itemsize != sizeof(float) ||
      pet_view.itemsize != sizeof(float) ||
      t1c_view.itemsize != sizeof(float) ||
      flair_view.itemsize != sizeof(float)) {
    PyErr_SetString(PyExc_ValueError, "MODEL, PET, T1c, FLAIR: wrong type");
    return NULL;
  }
  params.PET = PyBuffer_GetPointer(&pet_view, indices);
  params.T1c = PyBuffer_GetPointer(&t1c_view, indices);
  params.FLAIR = PyBuffer_GetPointer(&flair_view, indices);
  params.n[0] = model_view.shape[0];
  params.n[1] = model_view.shape[1];
  params.n[2] = model_view.shape[2];

  model = PyBuffer_GetPointer(&model_view, indices);
  params.PET = PyBuffer_GetPointer(&pet_view, indices);
  params.T1c = PyBuffer_GetPointer(&t1c_view, indices);
  params.FLAIR = PyBuffer_GetPointer(&flair_view, indices);
  brain_likelihood(&params, model, &ans);
  return PyFloat_FromDouble(ans);
}

static PyObject *run(PyObject *self, PyObject *args) {
  Py_buffer gm_view, wm_view, hg_view;
  Py_ssize_t indices[3] = {0, 0, 0};
  double dw, tend, h, t, L;
  int step, gpd;
  PyObject *GM, *WM, *HG;
  float *hg;
  struct Brain *brain;
  struct BrainParams params;

  if (!PyArg_ParseTuple(args, "iOO(ddd)dddO", &params.blocksPerDimension, &GM,
                        &WM, &params.ic[0], &params.ic[1], &params.ic[2], &dw,
                        &params.rho, &tend, &HG))
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
  // printf("%g %g %g\n", params.ic[0], params.ic[1], params.ic[2]);
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
  params.Dw = dw / (L * L);
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
  Py_RETURN_NONE;
}

static PyMethodDef Glioma_SolverMethods[] = {
    {"run", run, METH_VARARGS, PyDoc_STR("Run simulations")},
    {"likelihood", likelihood, METH_VARARGS, PyDoc_STR("Run simulations")},
    {NULL, NULL, 0, NULL},
};

static PyModuleDef glioma_solvermodule = {
    PyModuleDef_HEAD_INIT, "glioma_solver", PyDoc_STR("Glioma Solver"), -1,
    Glioma_SolverMethods,
};

PyMODINIT_FUNC PyInit_glioma_solver(void) {
  return PyModule_Create(&glioma_solvermodule);
}
