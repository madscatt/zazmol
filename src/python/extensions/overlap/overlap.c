#include <Python.h>
#include <math.h>

static PyObject* overlap(PyObject* self, PyObject* args) {
    PyObject *coor1, *coor2;
    double cut;
    int natoms1, natoms2, check = 0;
    if (!PyArg_ParseTuple(args, "OOid", &coor1, &coor2, &natoms1, &natoms2, &cut)) {
        return NULL;
    }

    for (int i = 0; i < natoms1; i++) {
        double x1 = PyFloat_AsDouble(PyList_GetItem(coor1, i*3));
        double y1 = PyFloat_AsDouble(PyList_GetItem(coor1, i*3 + 1));
        double z1 = PyFloat_AsDouble(PyList_GetItem(coor1, i*3 + 2));

        for (int j = 0; j < natoms2; j++) {
            double x2 = PyFloat_AsDouble(PyList_GetItem(coor2, j*3));
            double y2 = PyFloat_AsDouble(PyList_GetItem(coor2, j*3 + 1));
            double z2 = PyFloat_AsDouble(PyList_GetItem(coor2, j*3 + 2));

            double diff2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
            double dist = sqrt(diff2);

            if (dist < cut) {
                check = 1;
                break;
            }
        }
        if (check == 1) {
            break;
        }
    }

    return Py_BuildValue("i", check);
}

static PyMethodDef OverlapMethods[] = {
    {"overlap", overlap, METH_VARARGS, "Check overlap between atoms."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef overlapmodule = {
    PyModuleDef_HEAD_INIT,
    "overlap",
    NULL,
    -1,
    OverlapMethods
};

PyMODINIT_FUNC PyInit_overlap(void) {
    return PyModule_Create(&overlapmodule);
}
