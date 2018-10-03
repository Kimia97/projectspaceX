#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 20:26:05 2018

@author: rivertz
"""

import numpy as np
import math as m
import sys
import time
import matplotlib.pyplot as plot

lokalfeil = 0

class RungeKuttaFehlberg54:
    lokalfeil =0

    A = np.array(
        [[0, 0, 0, 0, 0, 0],
         [1 / 4, 0, 0, 0, 0, 0],
         [3 / 32, 9 / 32, 0, 0, 0, 0],
         [1932 / 2197, -7200 / 2197, 7296 / 2197, 0, 0, 0],
         [439 / 216, -8, 3680 / 513, -845 / 4104, 0, 0],
         [-8 / 27, 2, -3544 / 2565, 1859 / 4104, -11 / 40, 0]])

    B = np.array(
        [[25 / 216, 0, 1408 / 2565, 2197 / 4104, -1 / 5, 0],
         [16 / 135, 0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55]])

    def __init__(self,
                 function,
                 dimension,
                 stepsize,
                 tolerance):
        self.F = function
        self.dim = dimension
        self.h = stepsize
        self.tol = tolerance
        self.original_h = stepsize

    def step(self, Win):
        s = np.zeros((6, self.dim))

        for i in range(0, 6):
            s[i, :] = self.F(Win + self.h * self.A[i, 0:i].dot(s[0:i, :])) #Hva skjer her?

        Zout = Win + self.h * (self.B[0, :].dot(s)) #Sjekk steg for orden 5
        Wout = Win + self.h * (self.B[1, :].dot(s))

        E = np.linalg.norm(Wout - Zout, 2) / np.linalg.norm(Wout, 2) #Feilestimat
        return Wout, E

    def safeStep(self, Win):
        Wout, E = self.step(Win)
        # Check if the error is tolerable

        if (not self.isErrorTolerated(E)): #Dersom feilen er høyere enn toleransen fra bruker
            self.adjustStep(E)  #Endre trinnstørrelse for å bli optimal
            Wout, E = self.step(Win) #gjenta steg
        # hvis steget fortsatt ikke er innenfor toleransen
        counter = 0

        while not self.isErrorTolerated(E): #Når feilen er over toleransen
            self.divideStepByTwo() #Dele steg på to
            Wout, E = self.step(Win)
            counter = counter + 1
            if counter > 10:
                sys.exit(-1)

        self.adjustStep(E) #Endre stegstørrelse igjen
        self.lokalfeil += E
        return Wout, E

    def isErrorTolerated(self, E):
        return E < self.tol

    def adjustStep(self, E):
        if (E == 0):
            s = 2
        else:
            s = m.pow(self.tol * self.h / (2 * E), 0.25) #Pass på at E ikke blir null, kan bli feil her

        self.h = s * self.h
        if self.h > self.original_h:
            self.h = self.original_h


    def divideStepByTwo(self):
        self.h = self.h / 2

    def setStepLength(self, stepLength):
        self.h = stepLength



# Her skriver du inn initialverdiproblemet
def F(Y):
    res = np.zeros(3)
    res[0] = 1
    res[1] = -Y[1] - Y[2]
    res[2] = Y[1] - Y[2]
    return res

def finn_best_tol(F, dim, h, tEnd, W, ant_tester):
    tol = 05e-14
    modber = np.zeros(ant_tester)
    toleranser = np.zeros(ant_tester)

    for i in range(0, ant_tester):
        W = np.array([0, 1, 0])
        rkf45 = RungeKuttaFehlberg54(F, dim, h, tol)
        start_time = time.time()
        while W[0] < tEnd:
            W, E = rkf45.safeStep(W)
        rkf45.setStepLength(tEnd - W[0])
        W, E = rkf45.step(W)
        stop_time = time.time()

        #print('start:', start_time)
        #print('stop:', stop_time)

        bertid = stop_time - start_time
        print("Tid" , bertid)
        modber[i] = tEnd/bertid
        # print(i, "Toleranse: ", tol, " Tid: ", tid, " sekunder")
        toleranser[i] = tol
        tol = tol*10

    plot.loglog(toleranser, modber)
    plot.xlabel("Toleranse")
    plot.ylabel("Modeltid/beregningstid")
    plot.title("Funksjon av toleranse")
    plot.grid(True)
    plot.show()

def global_feil(w_1, y_1):
    return abs(w_1 - y_1)

def main():
    W = np.array([0, 1, 0])
    h = 0.25
    tol = 0.5e-14
    tEnd = 1.0
    rkf54 = RungeKuttaFehlberg54(F, 3, h, tol)

    while W[0] < tEnd:
        W, E = rkf54.safeStep(W)

    rkf54.setStepLength(tEnd - W[0])
    W, E = rkf54.step(W)

    print(W, E)
    y_1 = (m.e**(-1))*m.cos(1)
    y_2 = (m.e**(-1))*m.sin(1)
    print("Y_1: ", y_1)
    print("W_1: ", W[1])
    print("Y_2: ", y_2)
    print("W_2 ", W[2])

    globalfeil = global_feil(W[1], y_1)
    globalfeil2 = global_feil(W[2], y_2)
    print("Globalfeil1: ", globalfeil)
    print("Globalfeil2: ", globalfeil2)
    print("Lokalfeil: ", rkf54.lokalfeil)

    ant_tester = 10
    finn_best_tol(F, 3, 1/30, 10000, W, ant_tester)


if __name__ == "__main__":
    # execute only if run as a script
    main()
