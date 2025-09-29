import os
import pathlib
import ctypes
import numpy as np
from pymoo.indicators.igd import IGD
import matplotlib.pyplot as plt

import mowao_py

MOWAO_FILE = 'mowao.so'
mowao = mowao_py.MOWAO_PY(pathlib.Path().absolute() / MOWAO_FILE)


def ZDT1(x, x_len, f):
    f1 = x[0]
    g = 1 + 9.0 / (x_len - 1) * np.sum(x[1:x_len])
    f2 = g * (1 - np.power((f1 / g), 0.5))
    f[0] = f1
    f[1] = f2


def run(mw):
    pf = np.loadtxt('ZDT1.pf')
    igd = IGD(pf)
    pop = mw[0].pop
    rep = mw[0].rep

    igd_record = []
    for i in range(mw[0].maxiter):
        mowao.population_update(mw)
        mowao.repository_update(mw, mw[0].pop.pop)

        objs = np.array([])
        for j in range(rep.end):
            obj = []
            for k in range(mw[0].nobj):
                obj.append(rep.repo[j].pl[k])
            objs = np.append(objs, obj)

        objs = objs.reshape((int(objs.size / mw[0].nobj), mw[0].nobj))

        igd_record.append(igd(objs))
        print(
            "Iteration:", i + 1,
            "Population:", pop.end,
            "Repository:", rep.end,
            "IGD:", igd_record[-1],
        )

    print(
        'best:', np.min(igd_record),
        '\nworst:', np.max(igd_record),
        '\nmean:', np.mean(igd_record),
    )

    plt.plot(objs[:, 0], objs[:, 1], 'o')
    plt.xlabel('objective 1')
    plt.ylabel('objective 2')
    plt.show()


def main():
    mw = mowao.mowao_new()

    mw[0].maxiter = 100
    mw[0].hpop = 100
    mw[0].nrepo = 80
    mw[0].bond_radius = 0.001
    mw[0].punch = 0.5
    mw[0].evaporate = 0.5
    mw[0].coef = 2

    mw[0].nobj = 2
    mw[0].nvar = 10

    f = mowao_py.PROBLEM_FUNC(ZDT1)
    mw[0].f = f
    mowao.mowao_alloc(mw)

    for i in range(mw[0].nvar):
        mw[0].lb[i] = 0.0
        mw[0].ub[i] = 1.0
        mw[0].vlb[i] = -0.15
        mw[0].vub[i] = 0.15

    mowao.mowao_init(mw)
    # mowao.mowao_run(mw, 1)
    run(mw)
    mowao.mowao_clean(mw)
    mowao.mowao_free(mw)


if __name__ == '__main__':
    main()
