## IMPORT ##

import threading
import time
import math
from math import pow as p
from pynput.keyboard import Key, Listener
import numpy as np
import pygame
import object 
import random

## INITALISATION

running = True
pygame.init()

objects = []

## FUNCTIONS ##

def distance2D(a,b):
    dist = math.sqrt(p(a[0]-b[0],2) + p(a[1]-b[1],2))
    return dist

def mag(a):
    mag = math.sqrt(p(a[0],2)+p(a[1],2))
    return mag

def toScreenCoords(a, cenPos, scale, scr):
    a = (a - cenPos)/scale
    width = scr.get_width()
    height = scr.get_height()

    return (-a[0]+(width/2), a[1]+(height/2))

def toPhysCoords(a, cenPos, scale, scr):
    width = scr.get_width()
    height = scr.get_height()

    a -= np.array([width/2, height/2])
    a[0] = -a[0]
    a *= scale
    a += cenPos

    return a


def covalentForce(r, Eb, r0):
    r0 = r0/p(10,12)

    Eb = Eb/(6.02*p(10,26)) ##26

    k=Eb/p(r0,2)
    eTerm = np.exp(-k*(r-r0))
    modR = r/r0

    return 2*Eb*(1-eTerm/modR)*eTerm*(r*k+1)/(r0*p(modR,2))

def physicsLoop():
    global running
    global objects

    timeSpeed = 7*p(10,-11)
    ke = 8.9875 * p(10,9) 
    dt = 0.005

    while running:
        positions = []

        for i in objects:

            F = np.array([0,0])
            for j in objects:

                if i==j: 
                    pass
                else:
                    dist = distance2D(i.p,j.p)

                    ##COULOMBIC FORCE
                    Fc = np.array([0,0])
                    Fc = -(j.p-i.p)*ke*i.c*j.c/p(dist,3)

                    ##MORSE FORCES

                    Fmor = covalentForce(dist, 436, 62)*(j.p-i.p)/dist ##62

                    ##NET FORCE
                    F = F + Fmor

            a = F/i.m
            i.v = i.v + a*dt*timeSpeed
            positions.append(i.p + i.v*dt*timeSpeed)
        

        for i in range(len(positions)):
            
            objects[i].p = positions[i]
        
        time.sleep(dt)
    
    return

def displayLoop():
    global running
    global objects

    centerPosition = np.array([0,0])
    scale = p(10,-11) ##default is exponent -16

    scr = pygame.display.set_mode((600,500), pygame.RESIZABLE)
    pygame.display.set_caption('Moldyn')

    while running:
        pygame.display.flip()
        for event in pygame.event.get():  
            if event.type == pygame.QUIT:  
                running = False
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 4:
                    scale = scale/1.1
                elif event.button == 5:
                    scale = 1.1*scale
                elif event.button == 1:
                    mouse_x, mouse_y = pygame.mouse.get_pos()
                    mousePos = toPhysCoords(np.array([float(mouse_x), float(mouse_y)]), centerPosition, scale, scr)

                    objects.append(object.object(
                        5.3*p(10,-11),
                        mousePos,
                        0,
                        1.67*p(10, -27), 
                        0))


        keys = pygame.key.get_pressed()

        if keys[pygame.K_w]:
            centerPosition = centerPosition - scale*np.array([0,1])

        if keys[pygame.K_s]:
            centerPosition = centerPosition - scale*np.array([0,-1])

        if keys[pygame.K_d]:
            centerPosition = centerPosition - scale*np.array([1,0])

        if keys[pygame.K_a]:
            centerPosition = centerPosition - scale*np.array([-1,0])

        scr.fill((0,0,0))
        for i in objects:

            screenPosition = toScreenCoords(i.p, centerPosition, scale, scr)

            if i.c < 0 :
                pygame.draw.circle(scr, (0,0,255), screenPosition,i.r/scale)
            elif i.c > 0:
                pygame.draw.circle(scr, (255,0,0), screenPosition,i.r/scale)
            else:
                pygame.draw.circle(scr, (100,100,100), screenPosition,i.r/scale)
        
## CORE ##

if __name__ == '__main__':
    DThread = threading.Thread(target=displayLoop)
    PThread = threading.Thread(target=physicsLoop)

    DThread.start()
    PThread.start()
    DThread.join()
    PThread.join()

pygame.quit()
