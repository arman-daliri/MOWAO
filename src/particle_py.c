#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>

#include "types_py.h"
#include "../mowao/mowao.h"

PyObject *
particle_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	particle *self;

	self = (particle *)type->tp_alloc(type, 0);
	if (self == NULL)
		return NULL;

	self->dec = mwlist_py_type.tp_new(&mwlist_py_type, args, kwds);
	self->velocity = mwlist_py_type.tp_new(&mwlist_py_type, args, kwds);
	self->f = mwlist_py_type.tp_new(&mwlist_py_type, args, kwds);
	return (PyObject *)self;
}

void
particle_py_dealloc(PyObject *op)
{
	particle *self = (particle *)op;
	mwlist_py_dealloc(self->dec);
	mwlist_py_dealloc(self->velocity);
	mwlist_py_dealloc(self->f);
	Py_TYPE(self)->tp_free(self);
}

int
particle_py_module_exec(PyObject *m)
{
	if (PyType_Ready(&particle_py_type) < 0)
		return -1;
	if (PyModule_AddObjectRef(m, "particle",
				  (PyObject *)&particle_py_type) < 0)
		return -1;
	return 0;
}

static PyMemberDef particle_py_members[] = {
	{ "dec", Py_T_OBJECT_EX, offsetof(particle, dec), 0, "" },
	{ "velocity", Py_T_OBJECT_EX, offsetof(particle, velocity), 0, "" },
	{ "f", Py_T_OBJECT_EX, offsetof(particle, f), 0, "" },
	{ "evaporated", Py_T_INT, offsetof(particle, evaporated), 0, "" },
	{ "num_bond", Py_T_INT, offsetof(particle, num_bond), 0, "" },
	{ NULL },
};

// clang-format off
PyTypeObject particle_py_type = {
	.ob_base = PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "mowao.particle",
	.tp_doc = PyDoc_STR("particle type"),
	.tp_basicsize = sizeof(particle),
	.tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_new = particle_py_new,
  .tp_dealloc = particle_py_dealloc,
  .tp_members =particle_py_members,
};
// clang-format on
