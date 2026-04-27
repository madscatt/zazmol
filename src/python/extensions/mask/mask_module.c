#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

extern void get_mask_array(long long *farray, int nflexible, int natoms,
                           char **nname, long long *resid,
                           long long *flexible_residues, int nresidues,
                           int mtype);

static PyObject *
mask_get_mask_array(PyObject *self, PyObject *args)
{
    PyObject *farray_obj = NULL;
    PyObject *names_obj = NULL;
    PyObject *resid_obj = NULL;
    PyObject *flexible_obj = NULL;
    int nresidues = 0;
    int mtype = 0;

    PyArrayObject *farray = NULL;
    PyArrayObject *resid = NULL;
    PyArrayObject *flexible_residues = NULL;
    char **names = NULL;

    if (!PyArg_ParseTuple(args, "OOOOii", &farray_obj, &names_obj,
                          &resid_obj, &flexible_obj, &nresidues, &mtype)) {
        return NULL;
    }

    farray = (PyArrayObject *)PyArray_FROM_OTF(
        farray_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY);
    if (farray == NULL) {
        return NULL;
    }

    if (PyArray_NDIM(farray) != 2) {
        PyErr_Format(PyExc_TypeError,
                     "Array must have 2 dimensions.  Given array has %d dimensions",
                     PyArray_NDIM(farray));
        goto fail;
    }

    if (!PyArray_ISCARRAY(farray)) {
        PyErr_SetString(PyExc_TypeError,
                        "Array must be contiguous.  A non-contiguous array was given");
        goto fail;
    }

    if (PyArray_TYPE(farray) != NPY_LONGLONG) {
        PyErr_Format(PyExc_TypeError,
                     "Array of type 'long long' required.  Array of type '%s' given",
                     PyArray_DESCR(farray)->typeobj->tp_name);
        goto fail;
    }

    if (!PyList_Check(names_obj)) {
        PyErr_SetString(PyExc_TypeError, "not a list");
        goto fail;
    }

    Py_ssize_t name_count = PyList_Size(names_obj);
    names = (char **)malloc((name_count + 1) * sizeof(char *));
    if (names == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    for (Py_ssize_t i = 0; i < name_count; i++) {
        PyObject *item = PyList_GetItem(names_obj, i);
        if (!PyUnicode_Check(item)) {
            PyErr_SetString(PyExc_TypeError,
                            "list must contain strings: bozo");
            goto fail;
        }
        names[i] = (char *)PyUnicode_AsUTF8(item);
        if (names[i] == NULL) {
            goto fail;
        }
    }
    names[name_count] = NULL;

    resid = (PyArrayObject *)PyArray_FROM_OTF(
        resid_obj, NPY_LONGLONG, NPY_ARRAY_IN_ARRAY);
    if (resid == NULL) {
        goto fail;
    }

    flexible_residues = (PyArrayObject *)PyArray_FROM_OTF(
        flexible_obj, NPY_LONGLONG, NPY_ARRAY_IN_ARRAY);
    if (flexible_residues == NULL) {
        goto fail;
    }

    get_mask_array((long long *)PyArray_DATA(farray),
                   (int)PyArray_DIM(farray, 0),
                   (int)PyArray_DIM(farray, 1),
                   names,
                   (long long *)PyArray_DATA(resid),
                   (long long *)PyArray_DATA(flexible_residues),
                   nresidues,
                   mtype);

    free(names);
    Py_DECREF(flexible_residues);
    Py_DECREF(resid);
    PyArray_ResolveWritebackIfCopy(farray);
    Py_DECREF(farray);
    Py_RETURN_NONE;

fail:
    if (names != NULL) {
        free(names);
    }
    Py_XDECREF(flexible_residues);
    Py_XDECREF(resid);
    if (farray != NULL) {
        PyArray_DiscardWritebackIfCopy(farray);
        Py_DECREF(farray);
    }
    return NULL;
}

static PyMethodDef mask_methods[] = {
    {"get_mask_array", mask_get_mask_array, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mask_module = {
    PyModuleDef_HEAD_INIT,
    "_mask",
    NULL,
    -1,
    mask_methods
};

PyMODINIT_FUNC
PyInit__mask(void)
{
    PyObject *module = NULL;

    import_array();

    module = PyModule_Create(&mask_module);
    return module;
}
