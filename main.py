import threading
import time
import math
from math import pow as p
import numpy as np
import pygame
import object 
import random

pygame.init()

dt = 0.005
timeSpeed = 7*p(10,-11)
G = 6.67*p(10,-11) ##exponent -11
ke = 8.9875 * p(10,9) 
scale = 5*p(10,-12) ##default is exponent -16

objects = []
for i in range(16):
     objects.append(object.object(
        5.3*p(10,-11),
        p(10,-12)*np.array([random.randint(-100,100), random.randint(-100,100)]),
        0*p(10,-1)*np.array([random.randint(-10,30), random.randint(-10,30)]),
        1.67*p(10, -27), 
        0))

scr = pygame.display.set_mode((600,500), pygame.RESIZABLE)
pygame.display.set_caption('Moldyn')

def distance2D(a,b):
    dist = math.sqrt(p(a[0]-b[0],2) + p(a[1]-b[1],2))
    return dist

def mag(a):
    mag = math.sqrt(p(a[0],2)+p(a[1],2))
    return mag

def solvePDE():
    positions = []

    for i in objects:

        F = np.array([0,0])
        for j in objects:
            if i==j: pass
            else:
                dist = distance2D(i.p,j.p)

                ##COULOMBIC FORCE
                Fc = np.array([0,0])
                Fc = -(j.p-i.p)*ke*i.c*j.c/p(dist,3)

                ##MORSE FORCES

                Fmor = covalentForce(dist, 436, 62)*(j.p-i.p)/dist ##62

                ##NET FORCE
                F = F + Fmor
        print(mag(F))

        a = F/i.m
        i.v = i.v + a*dt*timeSpeed
        positions.append(i.p + i.v*dt*timeSpeed)
    
    for i in range(len(positions)):
        objects[i].p = positions[i]

def toScreenCoords(a):
    a = a/scale
    width = scr.get_width()
    height = scr.get_height()

    return (-a[0]+width/2, a[1]+height/2)

def covalentForce(r, Eb, r0):
    r0 = r0/p(10,12)

    Eb = Eb/(6.02*p(10,26)) ##26

    k=Eb/p(r0,2)
    eTerm = np.exp(-k*(r-r0))
    modR = r/r0

    return 2*Eb*(1-eTerm/modR)*eTerm*(r*k+1)/(r0*p(modR,2))

def rearrange(arr):
    for i in arr:
        j = random.randint(0,len(arr)-1)
        i = arr.index(i)
        e1 = arr[i]
        e2 = arr[j]
        arr[i] = e2
        arr[j] = e1

    return arr
        

running = True 
while running:

    time.sleep(dt)
    
    solvePDE()
    
    pygame.display.flip()
    for event in pygame.event.get():  
        if event.type == pygame.QUIT:  
            running = False

    scr.fill((0,0,0))
    for i in objects:

        screenPosition = toScreenCoords(i.p)

        if i.c < 0 :
            pygame.draw.circle(scr, (0,0,255), screenPosition,i.r/scale)
        elif i.c > 0:
            pygame.draw.circle(scr, (255,0,0), screenPosition,i.r/scale)
        else:
            pygame.draw.circle(scr, (100,100,100), screenPosition,i.r/scale)
        

pygame.quit()