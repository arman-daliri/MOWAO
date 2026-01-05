import numpy as np
from pymoo.indicators.igd import IGD
import matplotlib.pyplot as plt

import mowao


def ZDT1(x, f, cb_arg):
    f1 = x[0]
    g = 1 + 9.0 / (len(x) - 1) * sum(x[1:])
    f2 = g * (1 - pow((f1 / g), 0.5))
    f[0] = f1
    f[1] = f2


def run(mw):
    pf = np.loadtxt('ZDT1.pf')
    igd = IGD(pf)

    igd_record = []
    for i in range(mw.maxiter):
        mw.population_update()
        mw.repository_update()

        objs = np.array([])
        rep = mw.rep.list()
        for j in range(mw.rep.end):
            objs = np.append(objs, rep[j].f.list())

        objs = objs.reshape((int(objs.size / mw.nobj), mw.nobj))

        igd_record.append(igd(objs))
        print(
            "Iteration:", i + 1,
            "Population:", mw.pop.end,
            "Repository:", mw.rep.end,
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
    mw = mowao.mowao()
    mw.maxiter = 100
    mw.npop = 100
    mw.nrepo = 80
    mw.bond_radius = 0.001
    mw.push = 0.5
    mw.evaporate = 0.5
    mw.coef = 2

    mw.nobj = 2
    mw.ndec = 10
    mw.func_set(ZDT1, None)
    mw.alloc()
    for i in range(mw.ndec):
        mw.lb[i] = 0.0
        mw.ub[i] = 1.0
        mw.vlb[i] = -0.15
        mw.vub[i] = 0.15

    mw.init()
    run(mw)
    mw.clean()


main()
