#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>

#include "types_py.h"
#include "../mowao/mowao.h"

static int
carray_to_pylist(double *arr, PyObject *list, int len)
{
	PyObject *item;

	for (len--; -1 < len; len--) {
		item = PyFloat_FromDouble(arr[len]);
		if (item == NULL) {
			return -1;
		}
		if (PyList_SetItem(list, len, item) < 0) {
			return -1;
		}
	}
	return 0;
}

static int
pylist_to_carray(double *arr, PyObject *list, int len)
{
	double fd;
	PyObject *res;

	for (len--; -1 < len; len--) {
		res = PyList_GetItem(list, len);
		if (res == NULL)
			return -1;

		fd = PyFloat_AsDouble(res);
		arr[len] = fd;
	}
	return 0;
}

static void
f(double *x, int x_len, double *f, void *cb_arg)
{
	PyObject *res, *xlist = NULL, *flist = NULL;
	mowao *mw = (mowao *)cb_arg;

	xlist = PyList_New(x_len);
	flist = PyList_New(mw->mw.nobj);

	if (carray_to_pylist(x, xlist, x_len) < 0)
		goto cleanup;
	if (carray_to_pylist(f, flist, mw->mw.nobj) < 0)
		goto cleanup;

	res = PyObject_CallFunctionObjArgs(mw->py_callback, xlist, flist,
					   mw->py_cb_arg, NULL);
	Py_XDECREF(res);

	if (pylist_to_carray(f, flist, mw->mw.nobj) < 0)
		goto cleanup;
cleanup:
	Py_XDECREF(xlist);
	Py_XDECREF(flist);
}

static PyObject *
mowao_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	mowao *self;
	type = &mowao_py_type;
	self = (mowao *)type->tp_alloc(type, 0);
	self->lb = mwlist_py_type.tp_new(&mwlist_py_type, args, kwds);
	self->ub = mwlist_py_type.tp_new(&mwlist_py_type, args, kwds);
	self->vlb = mwlist_py_type.tp_new(&mwlist_py_type, args, kwds);
	self->vub = mwlist_py_type.tp_new(&mwlist_py_type, args, kwds);
	self->rep = repository_py_type.tp_new(&repository_py_type, args, kwds);
	self->pop = population_py_type.tp_new(&population_py_type, args, kwds);

	self->mw.cb_arg = (void *)self;
	self->mw.f = f;
	return (PyObject *)self;
}

static void
mowao_py_dealloc(PyObject *op)
{
	mowao *self = (mowao *)op;

	if (self->init) {
		population_clean(&(self->mw));
		repository_clean(&(self->mw));
	}
	mowao_clean_memebers(&(self->mw));
	if (self->alloc)
		mowao_clean_bounds(&(self->mw));
	Py_CLEAR(self->py_callback);
	Py_TYPE(self)->tp_free(self);
}

static int
mowao_py_init(PyObject *op, PyObject *args, PyObject *kwds)
{
	mowao *self = (mowao *)op;
	mwlist_py_type.tp_init(self->lb, args, kwds);
	mwlist_py_type.tp_init(self->ub, args, kwds);
	mwlist_py_type.tp_init(self->vlb, args, kwds);
	mwlist_py_type.tp_init(self->vub, args, kwds);
	repository_py_type.tp_init(self->rep, args, kwds);
	population_py_type.tp_init(self->pop, args, kwds);

	self->mw.cb_arg = (void *)self;
	self->mw.f = f;
	self->init = 0;
	self->alloc = 0;
	return 0;
}

static PyObject *
mowao_py_alloc(PyObject *op, PyObject *args)
{
	UNUSED(args);

	mowao *self = (mowao *)op;
	mowao_alloc(&(self->mw));
	mwlist *lb = (mwlist *)self->lb;
	mwlist *ub = (mwlist *)self->ub;
	mwlist *vlb = (mwlist *)self->vlb;
	mwlist *vub = (mwlist *)self->vub;

	lb->list = self->mw.lb;
	ub->list = self->mw.ub;
	vlb->list = self->mw.vlb;
	vub->list = self->mw.vub;

	lb->len = ub->len = self->mw.ndec;
	vlb->len = vub->len = self->mw.ndec;
	self->alloc = 1;
	return Py_None;
}

static PyObject *
mowao_py_initm(PyObject *op, PyObject *args)
{
	UNUSED(args);

	mowao *self = (mowao *)op;
	mowao_init(&(self->mw));
	repository *rep = (repository *)self->rep;
	population *pop = (population *)self->pop;

	rep->mw = &(self->mw);
	rep->rep = &(self->mw.rep);
	pop->mw = &(self->mw);
	pop->pop = &(self->mw.pop);

	self->init = 1;
	return Py_None;
}

static PyObject *
mowao_py_run(PyObject *op, PyObject *args)
{
	mowao *self = (mowao *)op;
	int log = 1;
	if (!PyArg_ParseTuple(args, "i", &log))
		return NULL;
	mowao_run(&(self->mw), log);
	return Py_None;
}

static PyObject *
mowao_py_clean(PyObject *op, PyObject *args)
{
	UNUSED(args);

	mowao *self = (mowao *)op;
	mowao_clean(&(self->mw));
	return Py_None;
}

static PyObject *
mowao_py_func_set(PyObject *op, PyObject *args)
{
	mowao *self = (mowao *)op;
	PyObject *temp;
	PyObject *py_cb_arg;

	if (PyArg_ParseTuple(args, "OO:set_callback;", &temp, &py_cb_arg)) {
		if (!PyCallable_Check(temp)) {
			PyErr_SetString(PyExc_TypeError,
					"parameter must be callable");
			return NULL;
		}
		Py_XINCREF(temp);
		Py_XDECREF(self->py_callback);
		self->py_callback = temp;
		self->py_cb_arg = py_cb_arg;
		return Py_None;
	}
	return NULL;
}

static PyObject *
mowao_py_population_update(PyObject *op, PyObject *args)
{
	UNUSED(args);

	mowao *self = (mowao *)op;
	population_update(&(self->mw));
	return Py_None;
}

static PyObject *
mowao_py_repository_update(PyObject *op, PyObject *args)
{
	UNUSED(args);

	mowao *self = (mowao *)op;
	repository_update(&(self->mw), self->mw.pop.pop);
	return Py_None;
}

int
mowao_py_module_exec(PyObject *m)
{
	if (PyType_Ready(&mowao_py_type) < 0)
		return -1;
	if (PyModule_AddObjectRef(m, "mowao", (PyObject *)&mowao_py_type) < 0)
		return -1;
	return 0;
}

static PyMethodDef mowao_py_methods[] = {
	{ "alloc", mowao_py_alloc, METH_NOARGS, "" },
	{ "init", mowao_py_initm, METH_NOARGS, "" },
	{ "run", mowao_py_run, METH_VARARGS, "" },
	{ "clean", mowao_py_clean, METH_NOARGS, "" },
	{ "func_set", mowao_py_func_set, METH_VARARGS, "" },
	{ "population_update", mowao_py_population_update, METH_NOARGS, "" },
	{ "repository_update", mowao_py_repository_update, METH_NOARGS, "" },
	{ NULL },
};

static PyMemberDef mowao_py_members[] = {
	{ "nobj", Py_T_INT, offsetof(mowao, mw.nobj), 0, "" },
	{ "ndec", Py_T_INT, offsetof(mowao, mw.ndec), 0, "" },
	{ "nrepo", Py_T_INT, offsetof(mowao, mw.nrepo), 0, "" },
	{ "npop", Py_T_INT, offsetof(mowao, mw.npop), 0, "" },
	{ "maxiter", Py_T_INT, offsetof(mowao, mw.maxiter), 0, "" },
	{ "bond_radius", Py_T_DOUBLE, offsetof(mowao, mw.bond_radius), 0, "" },
	{ "push", Py_T_DOUBLE, offsetof(mowao, mw.push), 0, "" },
	{ "evaporate", Py_T_DOUBLE, offsetof(mowao, mw.evaporate), 0, "" },
	{ "coef", Py_T_DOUBLE, offsetof(mowao, mw.coef), 0, "" },
	{ "lb", Py_T_OBJECT_EX, offsetof(mowao, lb), 0, "" },
	{ "ub", Py_T_OBJECT_EX, offsetof(mowao, ub), 0, "" },
	{ "vlb", Py_T_OBJECT_EX, offsetof(mowao, vlb), 0, "" },
	{ "vub", Py_T_OBJECT_EX, offsetof(mowao, vub), 0, "" },
	{ "rep", Py_T_OBJECT_EX, offsetof(mowao, rep), 0, "" },
	{ "pop", Py_T_OBJECT_EX, offsetof(mowao, pop), 0, "" },
	{ NULL },
};

// clang-format off
PyTypeObject mowao_py_type = {
	.ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "mowao.mowao",
	.tp_doc = PyDoc_STR("mowao type"),
	.tp_basicsize = sizeof(mowao),
	.tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT , 
  .tp_init = mowao_py_init,
  .tp_new = mowao_py_new,
  .tp_dealloc = mowao_py_dealloc,
  .tp_members = mowao_py_members,
  .tp_methods = mowao_py_methods,
};
// clang-format on
