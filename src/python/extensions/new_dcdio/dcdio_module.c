#include <Python.h>
#include <numpy/arrayobject.h>
#include "dcdio.h"  
#include <stdio.h>

int dum_read_dcdheader(FILE *fd, int *N, int *NSET, int *ISTART, int *NSAVC, double *DELTA, int *NAMNF, int *reverseEndian, int *charmm) ;

static PyObject* py_open_dcd_file(PyObject* self, PyObject* args) {
    const char* filename;

    // Parse the input tuple
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

// Function to open a DCD file for reading and return a capsule containing the file pointer
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


#define DCD_BADFORMAT -6

int my_reverseFourByteWord(int input) {
    return ((input >> 24) & 0x000000FF) |
           ((input >> 8)  & 0x0000FF00) |
           ((input << 8)  & 0x00FF0000) |
           ((input << 24) & 0xFF000000);
}

// Function to read the header of a DCD file
static PyObject* py_read_dcdheader(PyObject* self, PyObject* args) {
    PyObject* py_fd;
    int N, NSET, ISTART, NSAVC, NAMNF, reverseEndian, charmm;
    double DELTA;

    // Parse the input tuple
    if (!PyArg_ParseTuple(args, "O", &py_fd)) {
        fprintf(stderr, "Failed to parse arguments\n");
        fflush(stderr);
        return NULL;
    }

    FILE* fd = (FILE*)PyCapsule_GetPointer(py_fd, "dcdio_module.FILE");
    if (!fd) {
        PyErr_SetString(PyExc_ValueError, "Invalid file pointer");
        fprintf(stderr, "Invalid file pointer\n");
        fflush(stderr);
        return NULL;
    }

    // Ensure the file pointer is at the beginning of the file
    fseek(fd, 0, SEEK_SET);

    /*
    int input_integer;
    size_t ret_val = fread(&input_integer, sizeof(int), 1, fd);

    if (ret_val != 1) {
        PyErr_SetString(PyExc_IOError, "Error reading first int from DCD file");
        return NULL;
    }
    fprintf(stderr, "read_dcdheader: input_integer = %d\n", input_integer);
    fflush(stderr);

    // Check magic number in file header and determine byte order
    if (input_integer != 84) {
        // Reverse the byte order
        input_integer = my_reverseFourByteWord(input_integer);
        //input_integer= *reverseFourByteWord(&input_integer);

        if (input_integer != 84) {
            fprintf(stderr, "Invalid magic number in file header: %d\n", input_integer);
            fflush(stderr);
            PyErr_SetString(PyExc_ValueError, "Invalid magic number in file header");
            return NULL;
        }
    }

    fprintf(stderr, "Valid DCD file header\n");
    fflush(stderr);

    // Ensure the file pointer is at the beginning of the file
    fseek(fd, 0, SEEK_SET);
    */

    // Print debugging information before calling read_dcdheader
    fprintf(stderr, "Calling read_dcdheader with file pointer: %p\n", fd);
    fflush(stderr);

    // Call the actual read_dcdheader function from dcdio.h
    int result = read_dcdheader(fd, &N, &NSET, &ISTART, &NSAVC, &DELTA, &NAMNF, &reverseEndian, &charmm);
    if (result != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to read DCD header");
        fprintf(stderr, "read_dcdheader failed with result: %d\n", result);
        fflush(stderr);
        return NULL;
    }

    // Print debugging information to stderr
    fprintf(stderr, "read_dcdheader result: %d\n", result);
    fprintf(stderr, "N: %d, NSET: %d, ISTART: %d, NSAVC: %d, NAMNF: %d, DELTA: %f, reverseEndian: %d, charmm: %d\n",
            N, NSET, ISTART, NSAVC, NAMNF, DELTA, reverseEndian, charmm);
    fflush(stderr);  // Flush the output buffer

    // Return the results as a tuple
    return Py_BuildValue("iOiiiiidii", result, py_fd, N, NSET, ISTART, NSAVC, NAMNF, DELTA, reverseEndian, charmm);
}

// Example implementation of read_dcdheader function
int dum_read_dcdheader(FILE *fd, int *N, int *NSET, int *ISTART, int *NSAVC, double *DELTA, int *NAMNF, int *reverseEndian, int *charmm) {
    // Example implementation
    // Read the header from the DCD file and populate the provided pointers
    // Return 0 on success, non-zero on failure

    // For demonstration purposes, let's assume the header is read successfully
    *N = 1000;
    *NSET = 200;
    *ISTART = 0;
    *NSAVC = 1;
    *NAMNF = 0;
    *DELTA = 0.02;
    *reverseEndian = 0;
    *charmm = 1;

    return 0; // Success
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
    {"open_dcd_file", py_open_dcd_file, METH_VARARGS, "Open DCD file and return file pointer capsule"},
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