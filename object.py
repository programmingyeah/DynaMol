import numpy as np
from math import pow as p

chargeConstant = 1.6*p(10,-19) ##-19 exponent

## ALL VALUES IN SI UNITS EXCEPT FOR BONDING CONSTANTS: 
# BOND LENGTH: pm
# BOND DISSOCIATION ENERGY (BDE) : kJ/mol 

idTable = [
    (
        (5.3 * p(10,-11), 3.7 * p(10,-11)),
        1.67 * p(10, -27),
        (  ##BDE: H-ID
            436,
            428,
        )
    ),
    (
        (6.6 * p(10,-11), 6.6 * p(10,-11)),
        2.66 * p(10,-26),
        (
            428,
            498,
        ),

    )
]

class EConfig:

    def distribute(self, e):
        Kmax = 2
        Lmax = 8
        Mmax = 18
        Nmax = 32
        Omax = 50
        Pmax = 72

        for i in range(1,e):
            if e > Kmax:
                self.K = Kmax
                e-=2

                if e > Lmax:
                    self.L = Lmax
                    e-=Lmax

                    if e > Mmax:
                        self.M = Mmax
                        e-=Mmax

                        if e > Nmax:
                            self.N = Nmax
                            e-=Nmax

                            if e > Omax:
                                self.O = Omax
                                e-=Omax

                                if e > Pmax:
                                    self.P = Pmax
                                    e-=Pmax
                                else:
                                    self.P = e
                                    return

                            else:
                                self.O = e
                                return

                        else:
                            self.N = e
                            return

                    else:
                        self.M = e
                        return

                else:
                    self.L = e
                    return

            else:
                self.K = e
                return

    def __init__(self, e):
        self.K = 0
        self.L = 0
        self.M = 0
        self.N = 0
        self.O = 0
        self.P = 0

        self.distribute(e)

    def outerShell(self):
        Kmax = 2
        Lmax = 8
        Mmax = 18
        Nmax = 32
        Omax = 50
        Pmax = 72

        if self.K != Kmax:
            if self.K == 0:
                return "K"
            return "K"

        if self.L != Lmax:
            if self.L == 0:
                return "K"
            return "L"

        if self.M != Mmax:
            if self.M == 0:
                return "L"
            return "M"

        if self.N != Nmax:
            if self.M == 0:
                return "M"
            return "N"

        if self.O != Omax:
            if self.M == 0:
                return "M"
            return "O"

        if self.P != Pmax:
            if self.M == 0:
                return "O"
            return "P"

        def excite():

            LToN = {"K":1,"L":2,"M":3,"N":4,"O":5,"P":6}
            NToL = {1:"K",2:"L",3:"M",4:"N",5:"O",6:"P"}

            outerShell = self.outerShell()
            self[outerShell] -= 1
            self[NToL[LToN[outerShell]+1]] += 1

class object:
    def __init__(self, id, position, velocity):
        self.id = id
        id -= 1
        self.r = idTable[id][0]
        self.p = position
        self.m = idTable[id][1]
        self.v = velocity
        self.c = 0*chargeConstant
        self.BondCount = 0
        self.EConfig = EConfig(id+1)
