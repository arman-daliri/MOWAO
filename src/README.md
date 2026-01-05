# MOWAO 

Multi-objective WAO algorithm python package.

## Example
```python
import mowao
import numpy as np
import matplotlib.pyplot as plt


def ZDT1(x, f, cb):
    f1 = x[0]
    g = 1 + 9.0 / (len(x) - 1) * sum(x[1:])
    f2 = g * (1 - pow((f1 / g), 0.5))
    f[0] = f1
    f[1] = f2


def main():
    mw = mowao.mowao()
    mw.maxiter = 100
    mw.hpop = 100
    mw.nrepo = 80
    mw.bond_radius = 0.001
    mw.push = 0.5
    mw.evaporate = 0.5
    mw.coef = 2

    mw.nobj = 2
    mw.nvar = 10
    mw.func_set(ZDT1, None)
    mw.alloc()
    for i in range(mw.nvar):
        mw.lb[i] = 0.0
        mw.ub[i] = 1.0
        mw.vlb[i] = -0.15
        mw.vub[i] = 0.15

    mw.init()
    mw.run(0)  # mw.run(1)
    show(mw)
    mw.clean()


main()
```
