#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>

#include "types_py.h"
#include "../mowao/mowao.h"

PyObject *
population_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	UNUSED(args);
	UNUSED(kwds);

	population *self;
	self = (population *)type->tp_alloc(type, 0);
	if (self == NULL)
		return NULL;
	return (PyObject *)self;
}

void
population_py_dealloc(PyObject *op)
{
	population *self = (population *)op;
	Py_TYPE(self)->tp_free(self);
}

static PyObject *
population_py_getparticle(population *self, int ind)
{
	particle *p = (particle *)particle_py_type.tp_new(&particle_py_type,
							  NULL, NULL);
	mwlist *dec = (mwlist *)(p->dec);
	mwlist *velocity = (mwlist *)(p->velocity);
	mwlist *f = (mwlist *)(p->f);

	dec->list = self->pop->pop[ind].dec;
	velocity->list = self->pop->pop[ind].velocity;
	f->list = self->pop->pop[ind].f;

	dec->len = velocity->len = self->mw->ndec;
	f->len = self->mw->nobj;

	p->evaporated = self->pop->pop[ind].evaporated;
	p->num_bond = self->pop->pop[ind].num_bond;
	return (PyObject *)p;
}

static PyObject *
population_py_list(PyObject *op, PyObject *args)
{
	UNUSED(args);

	int i;
	PyObject *item;
	population *self = (population *)op;
	PyObject *l = PyList_New(self->pop->end);
	for (i = 0; i < self->pop->end; i++) {
		item = population_py_getparticle(self, i);
		Py_INCREF(item);
		PyList_SetItem(l, i, item);
	}
	return l;
}

static PyObject *
population_py_iter(PyObject *op)
{
	return PyObject_GetIter(population_py_list(op, NULL));
}

static PyObject *
population_py_end_get(PyObject *op, void *closure)
{
	UNUSED(closure);
	population *self = (population *)op;
	return PyLong_FromLong(self->pop->end);
}

static int
population_py_end_set(PyObject *op, PyObject *value, void *closure)
{
	UNUSED(op);
	UNUSED(value);
	UNUSED(closure);
	return -1;
}

int
population_py_module_exec(PyObject *m)
{
	if (PyType_Ready(&population_py_type) < 0)
		return -1;
	if (PyModule_AddObjectRef(m, "population",
				  (PyObject *)&population_py_type) < 0)
		return -1;
	return 0;
}

static PyMethodDef population_py_methods[] = {
	{ "list", population_py_list, METH_NOARGS, "" },
	{ NULL },
};

static PyGetSetDef population_py_getset[] = {
	{ "end", population_py_end_get, population_py_end_set, "", NULL },
	{ NULL },
};

// clang-format off
PyTypeObject population_py_type = {
	.ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "mowao.population",
	.tp_doc = PyDoc_STR("population type"),
	.tp_basicsize = sizeof(population),
	.tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_new = population_py_new,
  .tp_dealloc = population_py_dealloc,
  .tp_methods = population_py_methods,
  .tp_iter = population_py_iter,
  .tp_getset = population_py_getset,
};
// clang-format on
