#include <Python.h>
#include "dcdio.c"  // Include the dcdio.c file to access open_dcd_read

// Wrapper function for open_dcd_read
static PyObject* py_open_dcd_read(PyObject* self, PyObject* args) {
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    FILE* file = open_dcd_read(filename);
    if (file == NULL) {
        PyErr_SetString(PyExc_IOError, "Failed to open file");
        return NULL;
    }

    return PyLong_FromVoidPtr((void*)file);
}

// Method table
static PyMethodDef DcdMethods[] = {
    {"open_dcd_read", py_open_dcd_read, METH_VARARGS, "Open a DCD file for reading"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef dcdmodule = {
    PyModuleDef_HEAD_INIT,
    "dcdio_module",
    NULL,
    -1,
    DcdMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_dcdio_module(void) {
    return PyModule_Create(&dcdmodule);
}