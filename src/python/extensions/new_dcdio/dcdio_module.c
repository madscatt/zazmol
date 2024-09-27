#include <Python.h>
#include "dcdio.h"  // Assuming dcdio.h contains the declaration for read_dcdheader

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

// Module method definitions
static PyMethodDef DCDIOModuleMethods[] = {
    {"open_dcd_read", py_open_dcd_read, METH_VARARGS, "Open a DCD file for reading"},
    {"read_dcdheader", py_read_dcdheader, METH_VARARGS, "Read the header of a DCD file"},
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
    return PyModule_Create(&dcdio_module);
}