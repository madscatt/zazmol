#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include "dcdio.c"

// Wrapper for open_dcd_write
static PyObject* py_open_dcd_write(PyObject* self, PyObject* args) {
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }
    FILE* fd = open_dcd_write((char*)filename);  // Cast to char*
    if (fd == NULL) {
        PyErr_SetString(PyExc_IOError, "Failed to open file for writing");
        return NULL;
    }
    return PyCapsule_New(fd, "FILE*", NULL);
}

// Wrapper for close_dcd_write
static PyObject* py_close_dcd_write(PyObject* self, PyObject* args) {
    PyObject* capsule;
    if (!PyArg_ParseTuple(args, "O", &capsule)) {
        return NULL;
    }
    FILE* fd = (FILE*)PyCapsule_GetPointer(capsule, "FILE*");
    if (fd == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid file descriptor");
        return NULL;
    }
    int result = close_dcd_write(fd);
    return Py_BuildValue("i", result);
}

// Wrapper for pad
static PyObject* py_pad(PyObject* self, PyObject* args) {
    char* s;
    int len;
    if (!PyArg_ParseTuple(args, "si", &s, &len)) {
        return NULL;
    }
    pad(s, len);
    Py_RETURN_NONE;
}

// Wrapper for write_dcdheader
static PyObject* py_write_dcdheader(PyObject* self, PyObject* args) {
    PyObject* capsule;
    char* filename;
    int N, NSET, ISTART, NSAVC;
    double DELTA;
    if (!PyArg_ParseTuple(args, "Osiiiid", &capsule, &filename, &N, &NSET, &ISTART, &NSAVC, &DELTA)) {
        return NULL;
    }
    FILE* fd = (FILE*)PyCapsule_GetPointer(capsule, "FILE*");
    if (fd == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid file descriptor");
        return NULL;
    }
    int result = write_dcdheader(fd, filename, N, NSET, ISTART, NSAVC, DELTA);
    return Py_BuildValue("i", result);
}

// Wrapper for write_dcdstep
static PyObject* py_write_dcdstep(PyObject* self, PyObject* args) {
    PyObject* capsule;
    int N, curframe;
    PyArrayObject *X, *Y, *Z;
    if (!PyArg_ParseTuple(args, "OiiO!O!O!", &capsule, &N, &curframe, &PyArray_Type, &X, &PyArray_Type, &Y, &PyArray_Type, &Z)) {
        return NULL;
    }
    FILE* fd = (FILE*)PyCapsule_GetPointer(capsule, "FILE*");
    if (fd == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid file descriptor");
        return NULL;
    }
    if (PyArray_TYPE(X) != NPY_FLOAT32 || PyArray_TYPE(Y) != NPY_FLOAT32 || PyArray_TYPE(Z) != NPY_FLOAT32) {
        PyErr_SetString(PyExc_TypeError, "Arrays must be of type float32");
        return NULL;
    }
    float* dataX = (float*)PyArray_DATA(X);
    float* dataY = (float*)PyArray_DATA(Y);
    float* dataZ = (float*)PyArray_DATA(Z);
    int result = write_dcdstep(fd, N, dataX, dataY, dataZ, curframe);
    return Py_BuildValue("i", result);
}

// Wrapper for open_dcd_read
static PyObject* py_open_dcd_read(PyObject* self, PyObject* args) {
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }
    FILE* fd = open_dcd_read((char*)filename);  // Cast to char*
    if (fd == NULL) {
        PyErr_SetString(PyExc_IOError, "Failed to open file for reading");
        return NULL;
    }
    return PyCapsule_New(fd, "FILE*", NULL);
}

// Wrapper for reverseFourByteWord
static PyObject* py_reverseFourByteWord(PyObject* self, PyObject* args) {
    int value;
    if (!PyArg_ParseTuple(args, "i", &value)) {
        return NULL;
    }
    int* result_ptr = reverseFourByteWord(&value);  // Pass address of value
    int result = *result_ptr;  // Dereference the pointer to get the integer value
    return Py_BuildValue("i", result);
}

// Wrapper for reverseEightByteWord
static PyObject* py_reverseEightByteWord(PyObject* self, PyObject* args) {
    double value;
    if (!PyArg_ParseTuple(args, "d", &value)) {
        return NULL;
    }
    double* result_ptr = reverseEightByteWord(&value);  // Pass address of value
    double result = *result_ptr;  // Dereference the pointer to get the double value
    return Py_BuildValue("d", result);
}

// Wrapper for read_dcdheader
static PyObject* py_read_dcdheader(PyObject* self, PyObject* args) {
    PyObject* capsule;
    int N, NSET, ISTART, NSAVC, NAMNF;
    int reverseEndian = 0; // or appropriate value
    int charmm = 1; // or appropriate value

    double DELTA;
    PyArrayObject* array;
    if (!PyArg_ParseTuple(args, "OiiiiifO!", &capsule, &N, &NSET, &ISTART, &NSAVC, &NAMNF, &DELTA, &PyArray_Type, &array)) {
        return NULL;
    }
    FILE* fd = (FILE*)PyCapsule_GetPointer(capsule, "FILE*");
    if (fd == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid file descriptor");
        return NULL;
    }
    if (PyArray_TYPE(array) != NPY_FLOAT32) {
        PyErr_SetString(PyExc_TypeError, "Array must be of type float32");
        return NULL;
    }
    float* data = (float*)PyArray_DATA(array);

    // Add the missing argument (e.g., a placeholder value)
    int extra_arg = 0;  // Adjust this as necessary based on the actual function definition

    int result = read_dcdheader(fd, &N, &NSET, &ISTART, &NSAVC, &NAMNF, &DELTA, data, &extra_arg, &reverseEndian, &charmm);
    return Py_BuildValue("i", result);
}

// Wrapper for read_dcdstep
static PyObject* py_read_dcdstep(PyObject* self, PyObject* args) {
    PyObject* capsule;
    int N, curframe, arg1, arg2, arg3;
    PyArrayObject *X, *Y, *Z;
    if (!PyArg_ParseTuple(args, "OiO!O!O!iiii", &capsule, &N, &PyArray_Type, &X, &PyArray_Type, &Y, &PyArray_Type, &Z, &curframe, &arg1, &arg2, &arg3)) {
        return NULL;
    }
    FILE* fd = (FILE*)PyCapsule_GetPointer(capsule, "FILE*");
    if (fd == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid file descriptor");
        return NULL;
    }
    if (PyArray_TYPE(X) != NPY_FLOAT32 || PyArray_TYPE(Y) != NPY_FLOAT32 || PyArray_TYPE(Z) != NPY_FLOAT32) {
        PyErr_SetString(PyExc_TypeError, "Arrays must be of type float32");
        return NULL;
    }
    float* dataX = (float*)PyArray_DATA(X);
    float* dataY = (float*)PyArray_DATA(Y);
    float* dataZ = (float*)PyArray_DATA(Z);
    int result = read_dcdstep(fd, N, dataX, dataY, dataZ, curframe, arg1, arg2, arg3);
    return Py_BuildValue("i", result);
}

// Wrapper for close_dcd_read
static PyObject* py_close_dcd_read(PyObject* self, PyObject* args) {
    PyObject* capsule;
    if (!PyArg_ParseTuple(args, "O", &capsule)) {
        return NULL;
    }
    FILE* fd = (FILE*)PyCapsule_GetPointer(capsule, "FILE*");
    if (fd == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid file descriptor");
        return NULL;
    }
    int result = close_dcd_read(fd);
    return Py_BuildValue("i", result);
}

// Method table
static PyMethodDef DcdioMethods[] = {
    {"open_dcd_write", py_open_dcd_write, METH_VARARGS, "Open DCD file for writing"},
    {"close_dcd_write", py_close_dcd_write, METH_VARARGS, "Close DCD file after writing"},
    {"pad", py_pad, METH_VARARGS, "Pad function"},
    {"write_dcdheader", py_write_dcdheader, METH_VARARGS, "Write DCD header"},
    {"write_dcdstep", py_write_dcdstep, METH_VARARGS, "Write DCD step"},
    {"open_dcd_read", py_open_dcd_read, METH_VARARGS, "Open DCD file for reading"},
    {"reverseFourByteWord", py_reverseFourByteWord, METH_VARARGS, "Reverse four-byte word"},
    {"reverseEightByteWord", py_reverseEightByteWord, METH_VARARGS, "Reverse eight-byte word"},
    {"read_dcdheader", py_read_dcdheader, METH_VARARGS, "Read DCD header"},
    {"read_dcdstep", py_read_dcdstep, METH_VARARGS, "Read DCD step"},
    {"close_dcd_read", py_close_dcd_read, METH_VARARGS, "Close DCD file after reading"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef dcdio_module = {
    PyModuleDef_HEAD_INIT,
    "dcdio",
    NULL,
    -1,
    DcdioMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_dcdio(void) {
    import_array();  // Initialize NumPy
    return PyModule_Create(&dcdio_module);
}
