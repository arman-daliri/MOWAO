#ifndef __PY_TYPES__
#define __PY_TYPES__

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "../mowao/mowao.h"

#define UNUSED(x) (void)(x)

// clang-format off
typedef struct {
	PyObject_HEAD
	double *list;
	int len;
} mwlist;

typedef struct {
	PyObject_HEAD
	PyObject *dec;
	PyObject *velocity;
	PyObject *f;
	int evaporated;
	int num_bond;
} particle;

typedef struct {
	PyObject_HEAD
	struct MOWAO *mw;
	struct Repository *rep;
} repository;

typedef struct {
	PyObject_HEAD
	struct MOWAO *mw;
	struct Population *pop;
} population;

typedef struct {
	PyObject_HEAD
	PyObject *lb, *ub;
	PyObject *vlb, *vub;
	PyObject *py_callback;
	PyObject *pop;
	PyObject *rep;
  PyObject *py_cb_arg;
	struct MOWAO mw;
  int init, alloc;
} mowao;
// clang-format on

/* mowao */
extern PyTypeObject mowao_py_type;
int mowao_py_module_exec(PyObject *m);

/* mwlist */
extern PyTypeObject mwlist_py_type;
PyObject *mwlist_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
void mwlist_py_dealloc(PyObject *op);
int mwlist_py_module_exec(PyObject *m);
PyObject *mwlist_py_list(PyObject *op, PyObject *args);

/* particle */
extern PyTypeObject particle_py_type;
PyObject *particle_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int particle_py_module_exec(PyObject *m);

/* repository */
extern PyTypeObject repository_py_type;
PyObject *repository_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
void repository_py_dealloc(PyObject *op);
int repository_py_module_exec(PyObject *m);

/* population */
extern PyTypeObject population_py_type;
PyObject *population_py_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
void population_py_dealloc(PyObject *op);
int population_py_module_exec(PyObject *m);

#endif /* __PY_TYPES__ */
