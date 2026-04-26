#include <Python.h>
#include <numpy/arrayobject.h>

#include "view_vmd.h"

static PyObject *py_send_coordinates_to_vmd(PyObject *self, PyObject *args) {
    PyObject *x_obj;
    PyObject *y_obj;
    PyObject *z_obj;
    int port;
    int flag;
    PyArrayObject *x_array = NULL;
    PyArrayObject *y_array = NULL;
    PyArrayObject *z_array = NULL;
    npy_intp len_x;
    npy_intp len_y;
    npy_intp len_z;
    int result;

    if (!PyArg_ParseTuple(args, "OOOip", &x_obj, &y_obj, &z_obj, &port, &flag)) {
        return NULL;
    }

    x_array = (PyArrayObject *)PyArray_FROM_OTF(
        x_obj, NPY_NOTYPE, NPY_ARRAY_IN_ARRAY);
    y_array = (PyArrayObject *)PyArray_FROM_OTF(
        y_obj, NPY_NOTYPE, NPY_ARRAY_IN_ARRAY);
    z_array = (PyArrayObject *)PyArray_FROM_OTF(
        z_obj, NPY_NOTYPE, NPY_ARRAY_IN_ARRAY);

    if (x_array == NULL || y_array == NULL || z_array == NULL) {
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(z_array);
        return NULL;
    }

    if (PyArray_NDIM(x_array) != 1 || PyArray_NDIM(y_array) != 1 || PyArray_NDIM(z_array) != 1) {
        PyErr_SetString(PyExc_ValueError, "Coordinate arrays must be one-dimensional");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(z_array);
        return NULL;
    }

    if (PyArray_TYPE(x_array) != NPY_FLOAT32 ||
        PyArray_TYPE(y_array) != NPY_FLOAT32 ||
        PyArray_TYPE(z_array) != NPY_FLOAT32) {
        PyErr_SetString(PyExc_TypeError, "Coordinate arrays must have dtype float32");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(z_array);
        return NULL;
    }

    len_x = PyArray_DIM(x_array, 0);
    len_y = PyArray_DIM(y_array, 0);
    len_z = PyArray_DIM(z_array, 0);

    if (len_x != len_y || len_x != len_z) {
        PyErr_Format(PyExc_ValueError, "Arrays of lengths (%zd,%zd,%zd) given",
                     len_x, len_y, len_z);
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(z_array);
        return NULL;
    }

    result = send_coordinates_to_vmd(
        (int)len_x,
        (float *)PyArray_DATA(x_array),
        (float *)PyArray_DATA(y_array),
        (float *)PyArray_DATA(z_array),
        port,
        flag);

    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(z_array);

    return PyLong_FromLong((long)result);
}

static PyMethodDef ViewVmdMethods[] = {
    {"send_coordinates_to_vmd", py_send_coordinates_to_vmd, METH_VARARGS,
     "Send one frame of coordinates to VMD."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef viewvmdmodule = {
    PyModuleDef_HEAD_INIT,
    "_view_vmd",
    NULL,
    -1,
    ViewVmdMethods
};

PyMODINIT_FUNC PyInit__view_vmd(void) {
    import_array();
    return PyModule_Create(&viewvmdmodule);
}
