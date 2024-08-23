#include <Python.h>
#include <numpy/arrayobject.h>

static void matrix_multiply(double *a, double *b, double *c, int dim_a1, int dim_a2, int dim_b2);

static PyObject* matrix_math(PyObject* self, PyObject* args) {
    PyObject *a_obj, *b_obj, *c_obj;
    int dim_a1, dim_a2, dim_b2;

    if (!PyArg_ParseTuple(args, "OOOiii", &a_obj, &b_obj, &c_obj, &dim_a1, &dim_a2, &dim_b2)) {
        return NULL;
    }

    double *a = (double *)PyArray_DATA((PyArrayObject *)a_obj);
    double *b = (double *)PyArray_DATA((PyArrayObject *)b_obj);
    double *c = (double *)PyArray_DATA((PyArrayObject *)c_obj);

    matrix_multiply(a, b, c, dim_a1, dim_a2, dim_b2);

    Py_RETURN_NONE;
}

static PyMethodDef MatrixMathMethods[] = {
    {"matrix_math", matrix_math, METH_VARARGS, "Multiply two matrices"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef matrixmathmodule= {
    PyModuleDef_HEAD_INIT,
    "matrix_math",
    NULL,
    -1,
    MatrixMathMethods
};

PyMODINIT_FUNC PyInit_matrix_math(void) {
    import_array();  // Necessary for NumPy
    return PyModule_Create(&matrixmathmodule);
}

static void matrix_multiply(double *a, double *b, double *c, int dim_a1, int dim_a2, int dim_b2) {
    int i, j, k;
    for (i = 0; i < dim_a1; i++) {
        for (j = 0; j < dim_b2; j++) {
            c[i * dim_b2 + j] = 0.0;
        }
    }

    for (i = 0; i < dim_a1; i++) {
        for (j = 0; j < dim_b2; j++) {
            for (k = 0; k < dim_a2; k++) {
                c[i * dim_b2 + j] += a[i * dim_a2 + k] * b[k * dim_b2 + j];
            }
        }
    }
}