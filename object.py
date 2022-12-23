import numpy as np
from math import pow as p

chargeConstant = 1.6*p(10,-19) ##-19 exponent

radii = {
    
}

class object:
    def __init__(self, size, position, velocity, mass, charge):
        self.r = size
        self.p = position
        self.m = mass
        self.v = velocity
        self.c = charge*chargeConstant