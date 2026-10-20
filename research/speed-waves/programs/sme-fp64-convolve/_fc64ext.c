/* _fc64ext.c -- a thin CPython extension over fc64.c (av-1, wave 15), so the per-call overhead is that of a C
 * function call (buffer protocol, GIL released during the work), not of ctypes.
 *   corr(x, g, o0, y, path) : y[t] = sum_k g[k] xp[o0+t+k], t < len(y); float64 or complex128 (by itemsize of y)
 * Build: clang -O3 -mcpu=apple-m4 -shared -undefined dynamic_lookup -I$PYINC _fc64ext.c fc64.c \
 *        -framework Accelerate -o _fc64ext$EXT_SUFFIX
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

int fc_corr_d(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, int path);
int fc_corr_z(const double *x, long nx, const double *g, int L, long o0, long nout, double *y, int path);
int fc_vdsp_d(const double *x, long nx, const double *g, int L, long o0, long nout, double *y);

static PyObject *corr(PyObject *self, PyObject *const *args, Py_ssize_t nargs) {
  if (nargs != 5) { PyErr_SetString(PyExc_TypeError, "corr(x, g, o0, y, path)"); return NULL; }
  Py_buffer bx, bg, by;
  if (PyObject_GetBuffer(args[0], &bx, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT) < 0) return NULL;
  if (PyObject_GetBuffer(args[1], &bg, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT) < 0) { PyBuffer_Release(&bx); return NULL; }
  if (PyObject_GetBuffer(args[3], &by, PyBUF_C_CONTIGUOUS | PyBUF_WRITABLE | PyBUF_FORMAT) < 0) {
    PyBuffer_Release(&bx); PyBuffer_Release(&bg); return NULL;
  }
  long o0 = PyLong_AsLong(args[2]);
  int path = (int)PyLong_AsLong(args[4]);
  int r = 0;
  int cplx = by.itemsize == 16;
  if (bx.itemsize != by.itemsize || bg.itemsize != by.itemsize || (by.itemsize != 8 && by.itemsize != 16)) {
    PyErr_SetString(PyExc_TypeError, "x, g, y must all be float64 or all complex128");
    r = -2;
  } else {
    long nx = bx.len / bx.itemsize, L = bg.len / bg.itemsize, nout = by.len / by.itemsize;
    Py_BEGIN_ALLOW_THREADS
    if (cplx) r = fc_corr_z(bx.buf, nx, bg.buf, (int)L, o0, nout, by.buf, path);
    else if (path == 9) r = fc_vdsp_d(bx.buf, nx, bg.buf, (int)L, o0, nout, by.buf);
    else r = fc_corr_d(bx.buf, nx, bg.buf, (int)L, o0, nout, by.buf, path);
    Py_END_ALLOW_THREADS
  }
  PyBuffer_Release(&bx); PyBuffer_Release(&bg); PyBuffer_Release(&by);
  if (r == -2) return NULL;
  if (r) { PyErr_NoMemory(); return NULL; }
  Py_RETURN_NONE;
}

static PyMethodDef methods[] = {{"corr", (PyCFunction)(void (*)(void))corr, METH_FASTCALL, "fc64 correlation primitive"},
                                {NULL, NULL, 0, NULL}};
static struct PyModuleDef mod = {PyModuleDef_HEAD_INIT, "_fc64ext", NULL, -1, methods};
PyMODINIT_FUNC PyInit__fc64ext(void) { return PyModule_Create(&mod); }
