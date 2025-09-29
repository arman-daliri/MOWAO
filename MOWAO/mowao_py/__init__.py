import ctypes
import pathlib


PROBLEM_FUNC = ctypes.CFUNCTYPE(
    None,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_double)
)


class Particle(ctypes.Structure):
    _fields_ = [
        ("h", ctypes.POINTER(ctypes.c_double)),
        ("velocity", ctypes.POINTER(ctypes.c_double)),
        ("pl", ctypes.POINTER(ctypes.c_double)),
        ("evaporated", ctypes.c_int),
        ("numh", ctypes.c_int),
    ]


class Repository(ctypes.Structure):
    _fields_ = [
        ("repo", ctypes.POINTER(Particle)),
        ("dist", ctypes.POINTER(ctypes.c_double)),
        ("max", ctypes.POINTER(ctypes.c_double)),
        ("end", ctypes.c_int),
    ]


class Population(ctypes.Structure):
    _fields_ = [
        ("pop", ctypes.POINTER(Particle)),
        ("end", ctypes.c_int),
    ]


class MOWAO(ctypes.Structure):
    _fields_ = [
        ("nobj", ctypes.c_int),
        ("nvar", ctypes.c_int),
        ("lb", ctypes.POINTER(ctypes.c_double)),
        ("ub", ctypes.POINTER(ctypes.c_double)),
        ("nrepo", ctypes.c_int),
        ("hpop", ctypes.c_int),
        ("maxiter", ctypes.c_int),
        ("bond_radius", ctypes.c_double),
        ("punch", ctypes.c_double),
        ("evaporate", ctypes.c_double),
        ("coef", ctypes.c_double),
        ("vlb", ctypes.POINTER(ctypes.c_double)),
        ("vub", ctypes.POINTER(ctypes.c_double)),
        ("f", PROBLEM_FUNC),
        ("pop", Population),
        ("rep", Repository)
    ]


class MOWAO_PY:
    def __init__(self, mowao_so_path='mowao.so'):
        self.c_lib = ctypes.CDLL(mowao_so_path)

        self.randlim = self.c_lib.randlim
        self.randlim.argtype = [ctypes.c_double, ctypes.c_double]
        self.randlim.restype = ctypes.c_double

        self.dominates = self.c_lib.dominates
        self.dominates.argtype = [
            ctypes.POINTER(
                ctypes.c_double), ctypes.POINTER(
                    ctypes.c_double), ctypes.c_int]
        self.dominates.restype = ctypes.c_int

        self.check_boundary = self.c_lib.check_boundary
        self.check_boundary.argtype = [
            ctypes.POINTER(
                ctypes.c_double), ctypes.c_int, ctypes.POINTER(
                    ctypes.c_double), ctypes.POINTER(
                        ctypes.c_double)]

        self.distance = self.c_lib.distance
        self.distance.argtype = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int, ctypes.POINTER(ctypes.c_double)
        ]
        self.distance.restype = ctypes.c_double

        self.mowao_new = self.c_lib.mowao_new
        self.mowao_new.restype = ctypes.POINTER(MOWAO)

        self.mowao_alloc = self.c_lib.mowao_alloc
        self.mowao_alloc.argtype = [ctypes.POINTER(MOWAO)]

        self.mowao_init = self.c_lib.mowao_init
        self.mowao_init.argtype = [ctypes.POINTER(MOWAO)]

        self.mowao_run = self.c_lib.mowao_run
        self.mowao_run.argtype = [ctypes.POINTER(MOWAO), ctypes.c_int]

        self.mowao_clean = self.c_lib.mowao_clean
        self.mowao_clean.argtype = [ctypes.POINTER(MOWAO)]

        self.mowao_free = self.c_lib.mowao_free
        self.mowao_free.argtype = [ctypes.POINTER(MOWAO)]

        self.particle_alloc = self.c_lib.particle_alloc
        self.particle_alloc.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle)]

        self.particle_init = self.c_lib.particle_init
        self.particle_init.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle)]

        self.particle_clean = self.c_lib.particle_clean
        self.particle_clean.argtype = [ctypes.POINTER(Particle)]

        self.particle_copy = self.c_lib.particle_copy
        self.particle_copy.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle),
            ctypes.POINTER(Particle)
        ]

        self.particle_punch = self.c_lib.particle_punch
        self.particle_punch.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle)]

        self.particle_evaporate = self.c_lib.particle_evaporate
        self.particle_evaporate.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle)]

        self.population_init = self.c_lib.population_init
        self.population_init.argtype = [ctypes.POINTER(MOWAO)]

        self.population_clean = self.c_lib.population_clean
        self.population_clean.argtype = [ctypes.POINTER(MOWAO)]

        self.population_update = self.c_lib.population_update
        self.population_update.argtype = [ctypes.POINTER(MOWAO)]

        self.population_fill = self.c_lib.population_fill
        self.population_fill.argtype = [ctypes.POINTER(MOWAO)]

        self.repository_init = self.c_lib.repository_init
        self.repository_init.argtype = [ctypes.POINTER(MOWAO)]

        self.repository_clean = self.c_lib.repository_clean
        self.repository_clean.argtype = [ctypes.POINTER(MOWAO)]

        self.repository_check = self.c_lib.repository_check
        self.repository_check.argtype = [ctypes.POINTER(MOWAO)]

        self.repository_add = self.c_lib.repository_add
        self.repository_add.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle)]

        self.repository_update = self.c_lib.repository_update
        self.repository_update.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle)]

        self.repository_dist = self.c_lib.repository_dist
        self.repository_dist.argtype = [
            ctypes.POINTER(MOWAO),
            ctypes.POINTER(Particle)]
        self.repository_dist.restype = ctypes.c_double
