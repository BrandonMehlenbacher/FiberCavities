# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 09:52:00 2019

@author: bmehl
"""

#this program will be used to get any and all FFPC parameters that we so desire
#this will be an ever expanding program hopefully and will include a lot of stuff from 
import math
import numpy as np
import matplotlib.pyplot as plt


class FiberCavity:
    def __init__(self, ROC1,ROC2,Diameter1,Diameter2,Transmission1,Transmission2,cavityLength):
        self.roc1 = ROC1
        self.roc2 = ROC2
        self.diameter1 = Diameter1
        self.diameter2 = Diameter2
        self.transmission1 = Transmission1
        self.transmission2 = Transmission2
        self.cavityLength = cavityLength
    def get_ROC1(self):
        return self.roc1
    def get_ROC2(self):
        return self.roc2
    def get_Diameter1(self):
        return self.diameter1
    def get_Diameter2(self):
        return self.diameter2
    def get_Transmission1(self):
        return self.transmission1
    def get_Transmission2(self):
        return self.transmission2
    def get_cavityLength(self):
        return self.transmission2
    def idealFinesse(self):
        return (math.pi)*(1-self.get_transmission1()/10**6)/(self.get_transmission1()/10**6)
    def realFinesse(self,absorb,scatt):
        try:
            return (math.pi/(absorb+scatt+self.get_transmission1()+self.get_transmission2()))*10**6
        except absorb < .01 or scatt < .01:
            print('absorbance and scattering need to be in ppm')
    def cavityDecayRate(self):
        return (3*10**8)*math.pi/(2*self.idealFinesse()*self.get_cavityLength())
    def beam_radius_properties(self,wl):
        R1 = self.get_ROC1()
        R2= self.get_ROC2()
        d = self.get_cavityLength()
        if R1 == np.inf  or R2 == np.inf:
            if R2 == np.inf:
                R2 = R1
                z1 = 0
                z2 = d
                z0 = (-d*(R2+d))**0.5
            else:
                z1 = (-d*(R2+d))/(R2+R1+2*d) #position of mirror 1
                z2 = z1+d # position of mirror 2
                z0 = ((-d*(R1+d)*(R2+d)*(R2+R1+d))/(R2+R1+2*d)**2)**.5 # rayleigh range resonator
                if isinstance(z0,complex): # in case the values for R1 and R2 are incorrect
                    R1 = -R1 
                    R2 = -R2
                    z1 = (-d*(R2+d))/(R2+R1+2*d)
                    z2 = z1+d
                    z0 = ((-d*(R1+d)*(R2+d)*(R2+R1+d))/(R2+R1+2*d)**2)**.5
        w0 = ((wl*z0)/np.pi)**0.5 # beam radius at the "center"
        w1 = w0*(1+(z1/z0)**2)**0.5 # beam radius at mirror 1
        w2 = w0*(1+(z2/z0)**2)**0.5 # beam radius at mirror 2
        values_array = np.array([w0,w1,w2,z0,z1,z2])
        return values_array
    def minimum_mode_volume(self,wl):
        w0 = self.beam_radius_properties(wl)[0]
        return (np.pi/4)*(w0**2)*(self.get_cavity_length())