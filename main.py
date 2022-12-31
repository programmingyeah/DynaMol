## IMPORT ##

import threading
import time
import math
from math import pow as p
import numpy as np
import pygame
import object
idData = object.idTable
import random
import decimal
import copy

## INITALISATION

pause = False
running = True
scale = p(10,-11)
bondingConstant = -1
pygame.init()

objects = []
formerBonds = []
bonds = []

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

    Eb = Eb/(6.02*p(10,20)) ##20

    recR = 1/r
    s = r0/p(2,1/6)
    sTerm = p(recR * s,6)

    return 24*Eb*(2*recR*p(sTerm,2)-recR*sTerm)

def repulsiveForce(r, Eb, r0):
    Eb = Eb/(6.02*p(10,20)) ##20

    recR = 1/r
    s = r0/p(2,1/6)
    sTerm = p(recR * s,6)

    return 48*Eb*recR*p(sTerm,2)

def covalentPotential(i,j, Eb, r0):
    r = distance2D(i.p,j.p)

    Eb = Eb/(6.02*p(10,20)) ##26

    s = r0/p(2,1/6)
    sTerm = p(s/r,6)

    return 4*Eb*(p(sTerm,2)-sTerm)

def areBonded(i,j, Eb, r0):
    
    r = distance2D(i.p,j.p)
    if r == 0: r = r0/100
    if r > 5*r0: return False ##26

    Epot = covalentPotential(i,j, Eb, r0)
    Eb = Eb/(6.02*p(10,20)) ##26
    Ekin = 0.5*i.m*p(np.linalg.norm(i.v-j.v),2)

    return Epot + Ekin < -0.02*Eb

t=0
def physicsLoop():
    global running
    global objects
    global pause
    global bonds
    global idData
    global scale
    global t

    timeSpeed = 2*p(10,-13)
    ke = 8.9875 * p(10,9) 
    dt = 0.001

    while running:
        if pause == False:
            formerBonds = bonds
            bonds = []

            for i in objects:
                for j in objects:
                    if i != j:
                        if (i.EConfig.valence() > 0) and (j.EConfig.valence() > 0):
                            if areBonded(i,j, idData[i.id-1][2][j.id-1], idData[i.id-1][0][1] + idData[j.id-1][0][1]):
                                if objects != []:
                                    
                                    foundPreviousBond = False
                                    for x,y in formerBonds:
                                        if (x==objects.index(i)):
                                            if not (y==objects.index(j)):
                                                foundPreviousBond = True
                                                break
                                        elif (y==objects.index(i)):
                                            if not (x==objects.index(j)):
                                                foundPreviousBond = True
                                                break

                                    if foundPreviousBond: break
                                    
                                    bonds.append((objects.index(i),objects.index(j)))
                                    i.EConfig.BondCount += 1
                                    j.EConfig.BondCount += 1

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

                        ##MORSE FORCES

                        Fmor = np.array([0,0])
                        Fphi = np.array([0,0])

                        if (j.EConfig.valence() > 0):
                            BDE = idData[i.id-1][2][j.id-1]
                            r0 = idData[i.id-1][0][1] + idData[j.id-1][0][1]
                            Fmor = covalentForce(dist, BDE, r0)*(i.p-j.p)/dist ##436
                            
                        else:

                            isBonded = False
                            for x,y in bonds:
                                if ((objects[x]==i) and (objects[y]==j)) or ((objects[y]==i) and (objects[x]==j)): isBonded = True

                            if isBonded:
                                BDE = idData[i.id-1][2][j.id-1]
                                r0 = idData[i.id-1][0][1] + idData[j.id-1][0][1]
                                Fmor = covalentForce(dist, BDE,r0)*(i.p-j.p)/dist ##436
                            else:
                                aU = 0.8854*52.9/(p(i.id,0.23)+p(j.id,0.23))
                                x = decimal.Decimal(dist*p(10,12)/aU)
                                try:
                                    ex = decimal.Decimal.exp(x)
                                except decimal.Overflow:
                                    ex = p(10,decimal.Decimal(9000))
                                BDE = idData[i.id-1][2][j.id-1]/(6.02*p(10,20))
                                r0 = idData[i.id-1][0][1] + idData[j.id-1][0][1]
                                Fphi = -repulsiveForce(dist, idData[i.id-1][2][j.id-1],idData[i.id-1][0][1] + idData[j.id-1][0][1])*(j.p-i.p)/dist + 1.602*p(10,-19)*(j.p-i.p)*(-3.2*0.1818*p(ex,-3.2)/aU - 0.9432*0.5099*p(ex,-0.9432)/aU - 0.4028*0.2802*p(ex,-0.4028)/aU - 0.2016*0.02817*p(ex,-0.2016)/aU)/dist

                        ##NET FORCE
                        F = F + Fmor + Fphi

                a = F/i.m
                i.v = i.v + a*dt*timeSpeed
                t = t + dt*timeSpeed
                r=np.linalg.norm(i.p)
                if r>=(5000*p(10,-12)-i.r[0]):
                    sphereToCenter = (-i.p)/np.linalg.norm(i.p)
                    angle = math.pi/2 - math.acos(np.dot(sphereToCenter,i.v/np.linalg.norm(i.v))) 
                    radial = np.linalg.norm(i.v)*math.sin(angle)*sphereToCenter
                    i.v = i.v - 2*radial

                positions.append(i.p + i.v*dt*timeSpeed)
            

            totalKE = 0
            for i in range(len(positions)):
                if objects != []: 
                    objects[i].EConfig.BondCount = 0
                    objects[i].p = positions[i]
                    totalKE += 0.5*objects[i].m*np.linalg.norm(objects[i].v)

            k_B = 1.38 * p(10,-23)
            k_I = 8.314*1000
            if len(objects) != 0:
                T = (2/3) * ((totalKE)/len(objects))/k_B
                print(str(T) + " K")
                P = (len(objects))/(6.02*p(10,23))*T*k_I/(math.pi*p(10,-14))
                print(str(P*2*math.pi*p(10,-7)*p(10,12)) + " pN")
                print(str(P*1000000)+" uPa")

        time.sleep(dt)
    
    return

def displayLoop():
    global running
    global objects
    global pause
    global bonds
    global t

    centerPosition = np.array([0,0])
    global scale

    scr = pygame.display.set_mode((600,500), pygame.RESIZABLE)
    pygame.display.set_caption('Moldyn')

    copyBeginPosition = None
    copyFinalPosition = None
    toPaste = []

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

                    id = 1
                    objects.append(object.object(
                        id,
                        mousePos,
                        np.array([0,0])))
                elif event.button == 3:
                    mouse_x, mouse_y = pygame.mouse.get_pos()
                    mousePos = toPhysCoords(np.array([float(mouse_x), float(mouse_y)]), centerPosition, scale, scr)

                    id = 8
                    objects.append(object.object(
                        id,
                        mousePos,
                        np.array([0,0])))
                elif event.button == 2:
                    mouse_x, mouse_y = pygame.mouse.get_pos()
                    mousePos = toPhysCoords(np.array([float(mouse_x), float(mouse_y)]), centerPosition, scale, scr)

                    id = 6
                    objects.append(object.object(
                        id,
                        mousePos,
                        np.array([0,0])))
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    pause = not pause
                elif event.key == pygame.K_r:
                    objects = []
                    bonds = []
                elif event.key == pygame.K_c:
                    print("nugguhgnugh")
                    copyFinalPosition = None
                    toPaste = []
                    mouse_x, mouse_y = pygame.mouse.get_pos()
                    copyBeginPosition = toPhysCoords(np.array([float(mouse_x), float(mouse_y)]), centerPosition, scale, scr)
                elif event.key == pygame.K_v:
                    copyFinalPosition = None
                    copyBeginPosition = None

                    for i in toPaste:
                        mouse_x, mouse_y = pygame.mouse.get_pos()
                        i = copy.copy(i)
                        i.p = toPhysCoords(np.array([float(mouse_x), float(mouse_y)]), centerPosition, scale, scr) + i.p
                        objects.append(i)
            elif event.type == pygame.KEYUP:
                if event.key == pygame.K_c:
                    mouse_x, mouse_y = pygame.mouse.get_pos()
                    copyFinalPosition = toPhysCoords(np.array([float(mouse_x), float(mouse_y)]), centerPosition, scale, scr)

                    for i in objects:
                        if (i.p[0] >= copyFinalPosition[0] and i.p[0] <= copyBeginPosition[0]) or (i.p[0] <= copyFinalPosition[0] and i.p[0] >= copyBeginPosition[0]):
                            if (i.p[1] >= copyFinalPosition[1] and i.p[1] <= copyBeginPosition[1]) or (i.p[1] <= copyFinalPosition[1] and i.p[1] >= copyBeginPosition[1]):
                                copied = copy.copy(i)
                                toPaste.append(copied)
                                toPaste[toPaste.index(copied)].p = copyBeginPosition - i.p



        keys = pygame.key.get_pressed()

        if keys[pygame.K_w]: 
            centerPosition = centerPosition - scale*np.array([0,1])

        if keys[pygame.K_s]:
            centerPosition = centerPosition - scale*np.array([0,-1])

        if keys[pygame.K_d]:
            centerPosition = centerPosition - scale*np.array([1,0])

        if keys[pygame.K_a]:
            centerPosition = centerPosition - scale*np.array([-1,0])

        scr.fill((30,30,30))
        pygame.draw.circle(scr, (0,0,0), toScreenCoords(np.array([0,0]), centerPosition, scale, scr), 5000*p(10,-12)/scale)

        for i,j in bonds:
            smallest = objects[i].r[0]
            if objects[i].r[0] > objects[j].r[0]: smallest = objects[j].r[0]
            pygame.draw.line(scr, (250,250,250), toScreenCoords(objects[i].p, centerPosition, scale, scr), toScreenCoords(objects[j].p, centerPosition, scale, scr), int(smallest/(2*scale)))
        
        for i in objects:

            screenPosition = toScreenCoords(i.p, centerPosition, scale, scr)

            if i.id == 1:
                pygame.draw.circle(scr, (250,250,250), screenPosition,i.r[0]/scale)
            elif i.id == 8:
                pygame.draw.circle(scr, (250,0,0), screenPosition,i.r[0]/scale)
            elif i.id == 6:
                pygame.draw.circle(scr, (60,60,60), screenPosition,i.r[0]/scale)
        
        if not(copyBeginPosition is None):
            if copyFinalPosition is None:
                A, B = pygame.mouse.get_pos()
                screenBeginPosition = toScreenCoords(copyBeginPosition, centerPosition, scale, scr)
                C, D = screenBeginPosition[0], screenBeginPosition[1]

                if B >= D:
                    if A >= C:
                        rect = pygame.Rect(C, D, A-C, B-D)
                    else:
                        rect = pygame.Rect(A, D, C-A, B-D)
                else:
                    if A >= C:
                        rect = pygame.Rect(C, B, A-C, D-B)
                    else:
                        rect = pygame.Rect(A, B, C-A, D-B)

                pygame.draw.rect(scr, (255,255,255), rect, 2)
            

        font = pygame.font.Font('freesansbold.ttf', 28)
        text_surface = font.render("time: "+str(math.floor(100*t/p(10,-12))/100)+"ps", 1, (250,250,250))
        scr.blit(text_surface, (100, 40))

## CORE ##

if __name__ == '__main__':
    DThread = threading.Thread(target=displayLoop)
    PThread = threading.Thread(target=physicsLoop)

    PThread.start()
    DThread.start()
    PThread.join()
    DThread.join()

pygame.quit()
