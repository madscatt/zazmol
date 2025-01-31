#include <Python.h>
#include <numpy/arrayobject.h>

// Function to multiply two matrices
static PyObject* matrix_multiply(PyObject* self, PyObject* args) {
    PyObject *a_obj, *b_obj;
    int dim_a1, dim_a2, dim_b2;

    // Parse the input tuple
    if (!PyArg_ParseTuple(args, "OOiii", &a_obj, &b_obj, &dim_a1, &dim_a2, &dim_b2)) {
        return NULL;
    }

    // Get pointers to the data as C-types
    float *a = (float *)PyArray_DATA((PyArrayObject *)a_obj);
    float *b = (float *)PyArray_DATA((PyArrayObject *)b_obj);

    // Create the result array
    npy_intp dims[2] = {dim_a1, dim_b2};
    PyObject *c_obj = PyArray_SimpleNew(2, dims, NPY_FLOAT);
    float *c = (float *)PyArray_DATA((PyArrayObject *)c_obj);

    // Perform matrix multiplication
    int i, j, k;
    for (i = 0; i < dim_a1; i++) {
        for (j = 0; j < dim_b2; j++) {
            c[i * dim_b2 + j] = 0.0;
            for (k = 0; k < dim_a2; k++) {
//                printf("a[%d][%d] = %f, b[%d][%d] = %f\n", i, k, a[i * dim_a2 + k], k, j, b[k * dim_b2 + j]);
                c[i * dim_b2 + j] += a[i * dim_a2 + k] * b[k * dim_b2 + j];
            }
 //           printf("c[%d][%d] = %f\n", i, j, c[i * dim_b2 + j]);
        }
    }

    return c_obj;
}

// Method definitions
static PyMethodDef MatrixMathMethods[] = {
    {"matrix_multiply", matrix_multiply, METH_VARARGS, "Multiply two matrices"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef matrixmathmodule = {
    PyModuleDef_HEAD_INIT,
    "matrix_math",
    NULL,
    -1,
    MatrixMathMethods
};

// Module initialization
PyMODINIT_FUNC PyInit_matrix_math(void) {
    import_array();  // Necessary for NumPy
    return PyModule_Create(&matrixmathmodule);
}
