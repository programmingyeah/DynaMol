import numpy as np
from math import pow as p

chargeConstant = 1.6*p(10,-19) ##-19 exponent

## ALL VALUES IN SI UNITS EXCEPT FOR BONDING CONSTANTS: 
# BOND DISSOCIATION ENERGY (BDE) : kJ/mol 

idTable = [
    (
        (5.3 * p(10,-11), 3.7 * p(10,-11)),
        1.67 * p(10, -27),
        (  ##BDE: H-ID
            436, ##436
            0,
            0,
            0,
            0,
            0,
            0,
            428,
        ),
        ()
    ),
    (),
    (),
    (),
    (),
    (),
    (),
    (
        (6.6 * p(10,-11), 6.6 * p(10,-11)),
        2.66 * p(10,-26),
        (
            428,
            0,
            0,
            0,
            0,
            0,
            0,
            498,
        ),
    )
]

class EConfig:

    def distribute(self, e):

        for i in range(0,5):
            if e > self.maxes[i]:
                self.Shells[i] = self.maxes[i]
                e-=self.maxes[i]
            else:
                self.Shells[i] = e
                return

    def __init__(self, e):
        self.Shells = [0,0,0,0,0,0]
        self.maxes = [2,8,18,32,50,72]

        self.BondCount = 0

        self.distribute(e)

    def outerShell(self):
        Kmax = 2
        Lmax = 8
        Mmax = 18
        Nmax = 32
        Omax = 50
        Pmax = 72

        for i in range(0,5):
            if self.Shells[i] != self.maxes[i]:
                if self.Shells[i] == 0:
                    if i==0: return i

                    return i-1
                return i

    def excite(self):

        outerShell = self.outerShell()
        self.Shells[outerShell] -= 1
        self[outerShell+1] += 1

    def valence(self):
        outerShell = self.outerShell()
        return self.maxes[outerShell] - (self.Shells[outerShell] + self.BondCount)

class object:
    def __init__(self, id, position, velocity):
        self.id = id
        id -= 1
        self.r = idTable[id][0]
        self.p = position
        self.m = idTable[id][1]
        self.v = velocity
        self.c = 0*chargeConstant
        self.EConfig = EConfig(id+1)
        self.EConfig.distribute(id+1)
