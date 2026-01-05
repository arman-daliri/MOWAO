#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>

#include "types_py.h"
#include "../mowao/mowao.h"

PyObject *
mwlist_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	UNUSED(args);
	UNUSED(kwds);

	mwlist *self;
	self = (mwlist *)type->tp_alloc(type, 0);
	if (self == NULL)
		return NULL;
	return (PyObject *)self;
}

void
mwlist_py_dealloc(PyObject *op)
{
	mwlist *self = (mwlist *)op;
	Py_TYPE(self)->tp_free(self);
}

static PyObject *
mwlist_py_getitem(PyObject *op, PyObject *key)
{
	double ikey;
	if (!PyArg_Parse(key, "d", &ikey))
		return Py_None;
	mwlist *self = (mwlist *)op;
	return PyFloat_FromDouble(self->list[(int)ikey]);
}

int
mwlist_py_setitem(PyObject *op, PyObject *key, PyObject *v)
{
	double ikey, iv;
	if (!PyArg_Parse(key, "d", &ikey))
		return -1;
	if (!PyArg_Parse(v, "d", &iv))
		return -1;

	mwlist *self = (mwlist *)op;
	self->list[(int)ikey] = iv;
	return 0;
}

PyObject *
mwlist_py_list(PyObject *op, PyObject *args)
{
	UNUSED(args);

	int i;
	mwlist *self = (mwlist *)op;
	PyObject *l = PyList_New(self->len);
	for (i = 0; i < self->len; i++) {
		PyObject *item = PyFloat_FromDouble(self->list[i]);
		Py_INCREF(item);
		PyList_SetItem(l, i, item);
	}
	return l;
}

static PyObject *
mwlist_py_iter(PyObject *op)
{
	return PyObject_GetIter(mwlist_py_list(op, NULL));
}

int
mwlist_py_module_exec(PyObject *m)
{
	if (PyType_Ready(&mwlist_py_type) < 0)
		return -1;
	if (PyModule_AddObjectRef(m, "mwlist", (PyObject *)&mwlist_py_type) < 0)
		return -1;
	return 0;
}

static PyMethodDef mwlist_py_methods[] = {
	{ "list", mwlist_py_list, METH_NOARGS, "" },
	{ NULL },
};

static PyMappingMethods mwlist_py_map = {
	.mp_subscript = mwlist_py_getitem,
	.mp_ass_subscript = mwlist_py_setitem,
};

PyObject *
new_alloc(PyTypeObject *type, Py_ssize_t nitems)
{
	PyObject *obj;
	const size_t size = _PyObject_VAR_SIZE(type, nitems + 1);

	obj = (PyObject *)malloc(size);
	if (obj == NULL) {
		PyErr_NoMemory();
		return NULL;
	}
	/*
 * If we don't need to zero memory, we could use
 * PyObject_{New, NewVar} for this whole function.
 */
	memset(obj, 0, size);
	if (type->tp_itemsize == 0) {
		PyObject_Init(obj, type);
	} else {
		(void)PyObject_InitVar((PyVarObject *)obj, type, nitems);
	}
	return obj;
}

// clang-format off
PyTypeObject mwlist_py_type = {
	.ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "mowao.mwlist",
	.tp_doc = PyDoc_STR("mwlist type"),
	.tp_basicsize = sizeof(mwlist),
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  
  .tp_new = mwlist_py_new,
  .tp_dealloc = mwlist_py_dealloc,
  .tp_methods = mwlist_py_methods,
  .tp_as_mapping = &mwlist_py_map,
  .tp_iter = mwlist_py_iter,
};
// clang-format on
