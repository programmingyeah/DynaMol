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
