#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>

#include "types_py.h"
#include "../mowao/mowao.h"

PyObject *
repository_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	UNUSED(args);
	UNUSED(kwds);

	repository *self;
	self = (repository *)type->tp_alloc(type, 0);
	if (self == NULL)
		return NULL;
	return (PyObject *)self;
}

void
repository_py_dealloc(PyObject *op)
{
	repository *self = (repository *)op;
	Py_TYPE(self)->tp_free(self);
}

static PyObject *
repository_py_getparticle(repository *self, int ind)
{
	particle *p = (particle *)particle_py_type.tp_new(&particle_py_type,
							  NULL, NULL);
	mwlist *dec = (mwlist *)(p->dec);
	mwlist *velocity = (mwlist *)(p->velocity);
	mwlist *f = (mwlist *)(p->f);

	dec->list = self->rep->repo[ind].dec;
	velocity->list = self->rep->repo[ind].velocity;
	f->list = self->rep->repo[ind].f;

	dec->len = velocity->len = self->mw->ndec;
	f->len = self->mw->nobj;

	p->evaporated = self->rep->repo[ind].evaporated;
	p->num_bond = self->rep->repo[ind].num_bond;

	return (PyObject *)p;
}

static PyObject *
repository_py_list(PyObject *op, PyObject *args)
{
	UNUSED(args);

	int i;
	PyObject *item;
	repository *self = (repository *)op;
	PyObject *l = PyList_New(self->rep->end);
	for (i = 0; i < self->rep->end; i++) {
		item = repository_py_getparticle(self, i);
		Py_INCREF(item);
		PyList_SetItem(l, i, item);
	}
	return l;
}

static PyObject *
repository_py_iter(PyObject *op)
{
	return PyObject_GetIter(repository_py_list(op, NULL));
}

static PyObject *
repository_py_end_get(PyObject *op, void *closure)
{
	UNUSED(closure);
	repository *self = (repository *)op;
	return PyLong_FromLong(self->rep->end);
}

static int
repository_py_end_set(PyObject *op, PyObject *value, void *closure)
{
	UNUSED(op);
	UNUSED(value);
	UNUSED(closure);
	return -1;
}

int
repository_py_module_exec(PyObject *m)
{
	if (PyType_Ready(&repository_py_type) < 0)
		return -1;
	if (PyModule_AddObjectRef(m, "repository",
				  (PyObject *)&repository_py_type) < 0)
		return -1;
	return 0;
}

static PyMethodDef repository_py_methods[] = {
	{ "list", repository_py_list, METH_NOARGS, "" },
	{ NULL },
};

static PyGetSetDef repository_py_getset[] = {
	{ "end", repository_py_end_get, repository_py_end_set, "", NULL },
	{ NULL },
};

// clang-format off
PyTypeObject repository_py_type = {
	.ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "mowao.repository",
	.tp_doc = PyDoc_STR("repository type"),
	.tp_basicsize = sizeof(repository),
	.tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_new = repository_py_new,
  .tp_dealloc = repository_py_dealloc,
  .tp_methods = repository_py_methods,
  .tp_iter = repository_py_iter,
  .tp_getset = repository_py_getset,
};
// clang-format on
