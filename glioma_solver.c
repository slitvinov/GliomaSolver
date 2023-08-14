#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject *run(PyObject *self, PyObject *args) {
  double x, y;
  if (!PyArg_ParseTuple(args, "dd", &x, &y))
    return NULL;
  return PyFloat_FromDouble(0.0);
}

static PyMethodDef Glioma_SolverMethods[] = {
    {"run", run, METH_VARARGS, PyDoc_STR("Run simulations")},
    {NULL, NULL, 0, NULL},
};

static PyModuleDef glioma_solvermodule = {
    PyModuleDef_HEAD_INIT,
    "glioma_solver",
    PyDoc_STR("Glioma Solver"),
    -1,
    Glioma_SolverMethods,
};

PyMODINIT_FUNC PyInit_glioma_solver(void) {
  return PyModule_Create(&glioma_solvermodule);
}
