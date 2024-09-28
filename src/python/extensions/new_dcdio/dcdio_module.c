#include <Python.h>
#include <numpy/arrayobject.h>
#include "dcdio.h"  

static PyObject* py_open_dcd_read(PyObject* self, PyObject* args) {
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    FILE* fd = fopen(filename, "rb");
    if (!fd) {
        PyErr_SetString(PyExc_IOError, "Failed to open file");
        return NULL;
    }

    return PyCapsule_New(fd, "dcdio_module.FILE", NULL);
}

static PyObject* py_read_dcdheader(PyObject* self, PyObject* args) {
    PyObject* py_fd;
    int N, NSET, ISTART, NSAVC, NAMNF, reverseEndian, charmm, extra_arg;
    double DELTA;
    float data;

    // Parse the input tuple
    if (!PyArg_ParseTuple(args, "Oiiiiidfiii", &py_fd, &N, &NSET, &ISTART, &NSAVC, &NAMNF, &DELTA, &data, &extra_arg, &reverseEndian, &charmm)) {
        return NULL;
    }

    FILE* fd = (FILE*)PyCapsule_GetPointer(py_fd, "dcdio_module.FILE");
    if (!fd) {
        PyErr_SetString(PyExc_ValueError, "Invalid file pointer");
        return NULL;
    }

    // Call the actual read_dcdheader function from dcdio.h
    int result = read_dcdheader(fd, &N, &NSET, &ISTART, &NSAVC, &NAMNF, &DELTA, &data, &extra_arg, &reverseEndian, &charmm);
    if (result != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to read DCD header");
        return NULL;
    }

    // Return the results as a tuple
    return Py_BuildValue("iiiiidfiii", N, NSET, ISTART, NSAVC, NAMNF, DELTA, data, extra_arg, reverseEndian, charmm);
}


// Function to read a DCD step

static PyObject* py_read_dcdstep(PyObject *self, PyObject *args) {
    PyObject *py_fp;
    int natoms, num_fixed, first, reverseEndian, charmm;
    PyArrayObject *x_array, *y_array, *z_array;

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "OiiiiiO!O!O!", &py_fp, &natoms, &num_fixed, &first, &reverseEndian, &charmm,
                          &PyArray_Type, &x_array, &PyArray_Type, &y_array, &PyArray_Type, &z_array)) {
        return NULL;
    }

    // Ensure the file pointer is a PyCapsule and extract the FILE* pointer
    FILE* fp = (FILE*) PyCapsule_GetPointer(py_fp, "dcdio_module.FILE");
    if (!fp) {
        PyErr_SetString(PyExc_ValueError, "Invalid file pointer");
        return NULL;
    }

    // Check if the arrays are NULL
    if (x_array == NULL || y_array == NULL || z_array == NULL) {
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(z_array);
        return NULL;
    }

    // Get pointers to the data as C-types
    float *x = (float*) PyArray_DATA(x_array);
    float *y = (float*) PyArray_DATA(y_array);
    float *z = (float*) PyArray_DATA(z_array);

    // Call the existing read_dcdstep function
    int result = read_dcdstep(fp, natoms, x, y, z, num_fixed, first, reverseEndian, charmm);

    if (result != 0) {  
        PyErr_SetString(PyExc_RuntimeError, "Failed to read DCD step");
        return NULL;
    }
    // Print the result to the screen
    printf("Result: %d\n", result);
    fflush(stdout);  // Flush the output buffer

    // Clean up
    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(z_array);

    // Return the result
    return Py_BuildValue("i", result);
}

// Module method definitions
static PyMethodDef DCDIOModuleMethods[] = {
    {"open_dcd_read", py_open_dcd_read, METH_VARARGS, "Open a DCD file for reading"},
    {"read_dcdheader", py_read_dcdheader, METH_VARARGS, "Read the header of a DCD file"},
    {"read_dcdstep", py_read_dcdstep, METH_VARARGS, "Read a step from a DCD file"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef dcdio_module = {
    PyModuleDef_HEAD_INIT,
    "dcdio_module",
    NULL,
    -1,
    DCDIOModuleMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_dcdio_module(void) {
    import_array();  // Initialize the NumPy API
    return PyModule_Create(&dcdio_module);
}