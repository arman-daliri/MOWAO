#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>

#include "types_py.h"
#include "../mowao/mowao.h"

static PyModuleDef_Slot slots[] = {
	{ Py_mod_exec, mowao_py_module_exec },
	{ Py_mod_exec, mwlist_py_module_exec },
	{ Py_mod_exec, particle_py_module_exec },
	{ Py_mod_exec, repository_py_module_exec },
	{ Py_mod_exec, population_py_module_exec },
	{ Py_mod_multiple_interpreters,
	  Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED },
	{ 0, NULL },
};

static struct PyModuleDef module = {
	.m_base = PyModuleDef_HEAD_INIT,
	.m_name = "mowao",
	.m_doc = "",
	.m_size = 0,
	.m_slots = slots,
};

PyMODINIT_FUNC
PyInit_mowao(void)
{
	return PyModuleDef_Init(&module);
}
