import numpy as np
import matplotlib.pyplot as plt
import time
import argparse

def efficiency(mode_radius,core_radius,radius_of_curvature,wavelength):
    """
    Inputs:
    mode_radius: radius of the cavity mode inside of the cavity assumed to be in microns
    core_radius: simple def, core of the field, complicated def, mode field radius of the field (extends slightly past core of fiber) assumed to be in microns
    radius_of_curvature: radius of curvature of the fiber, we are assuming a symmetric cavity ie both ROCs are the same. assumed to be in microns
    wavelength: wavelenght of light we plan to optimally use. Assumed to be in microns

    Outputs:
    efficiency: The output is the efficiency of coupling light into the fiber

    """
    #return ((2*mode_radius*core_radius)/(mode_radius**2+core_radius**2))**2
    return 4/(((mode_radius/core_radius)+(core_radius/mode_radius))**2+((np.pi*1.54*mode_radius*core_radius)/(wavelength*radius_of_curvature)))
def mode_radius_fn(radius_curvature,length_cavity,wavelength):
    """
    Inputs:
    radius_curvature: The radius of curvature of the mirror, this function assumes microns
    length_cavity: The length of the cavity, this function assumes length is in microns
    
    Outputs:
    mode_radius: The smallest the spot size gets within the cavity

    """
    return ((wavelength/(2*np.pi))**0.5)*(length_cavity*(2*radius_of_curvature-length_cavity))**0.25

def rayleigh_range(mode_radius_value, wavelength):
    """
    Inputs:
    spot_size: Smallest size inside of the cavity, assumed to be in microns
    wavelength: wavelength of light that will be used, assumed to be in microns

    Outputs:
    rayleigh_range: The rayleigh range if the system in microns

    """
    return (np.pi*mode_radius_value**2)/(wavelength)

def mode_radius_mirror_fn(cavity_length,mode_radius_value,wavelength):
    """
    Inputs:
    cavity_length: length of the cavity in microns
    mode_radius: the minimum mode waist inside of the cavity in microns
    wavelength: wavelength of interest in microns

    Outputs:
    mode_radius on the the cavity mirror

    """
    cavity_length = cavity_length/2
    rayleigh_range_value = rayleigh_range(mode_radius_value,wavelength)
    return mode_radius_value*(1+(cavity_length/rayleigh_range_value)**2)**0.5

def clipping_losses(radius_mirror,mode_radius_value):
    """
    Inputs:
    radius_mirror: radius of the mirror, this is half the diameter of the mirror if you go back to the old fiber spreadsheets, in microns
    mode_radius_value: the size of the mode when it reaches the mirror, calculated from previous mode_radius_mirror, in microns
    
    """
    return np.exp(-2*((radius_mirror)**2)/(mode_radius_value**2))
def finesse_calculations(cavity_length,mirror_diameter,mode_radius_value):
    """
    Inputs:
    cavity_length: the length of the cavity in microns
    mirror_diameter: the diameter of the mirror in microns
    mode_radius_value: the size of the mode when it reaches the mirror, calculated from previous mode_radius_mirror, in microns

    Outputs:
    Finesse: Calculates the finesse of the cavity in terms of losses
    """
    clipping_loss = clipping_losses(mirror_diameter,mode_radius_value)*2
    total_losses = 160*10**-6 # tbh i was lazy and decided a constant starting finesse value, this includes absorption,scattering, and transmission losses
    return (2*np.pi)/(clipping_loss+total_losses)

def cross_sectional_area(mode_radius):
    return np.pi*mode_radius

def diffraction_limited_spot(radius_of_curvature,mirror_diameter,wavelength):
    return (wavelength*radius_of_curvature)/(np.pi*mirror_diameter)
wavevector = lambda x: (2*np.pi)/(x) #x is wavelength and is assumed to be in microns
def cooperativity(finesse,wavevector,spot_size):
    return (24*finesse)/(np.pi*(wavevector**2)*spot_size**2)
#Defining Command Line arguments

savefig = False #decides whether or not to save the file
parser = argparse.ArgumentParser(description = "Takes in ROC and Mirror Diameter")
parser.add_argument("-m","--mirror_diameter", help = "Input mirror diameter in microns", required = False,default = 100)
parser.add_argument("-r","--radius_of_curvature", help = "Input radius of curvature in microns", required = False,default = 1000)
parser.add_argument("-s","--save",help = "If passed, enables the save, otherwise figure is displayed", required = False,default = "n")
parser.add_argument("-w","--wavelength",help = 'Wavelength of light in cavity in nanometers',required = False, default = 780)
parser.add_argument("-c","--core_fiber",help = 'This is the mode field diameter of an optical fiber in microns',required = False,default = 3)
argument = parser.parse_args()

if argument.mirror_diameter:
    arg_mir_dia = int(argument.mirror_diameter)
if argument.radius_of_curvature:
    arg_ROC = int(argument.radius_of_curvature)
if argument.save:
    if argument.save == "y":
        savefig = True
    else:
        savefig = False
if argument.wavelength:
    wavelength_nm = int(argument.wavelength)
if argument.core_fiber:
    core_radius = int(argument.core_fiber)


mirror_diameter_list = [arg_mir_dia]
ROC_list = [arg_ROC]
wavelength = wavelength_nm*10**-3
for radius_of_curvature in ROC_list:
    for diameter_mirror in mirror_diameter_list:
        radius_mirror = diameter_mirror/2
        max_cavity_length = radius_of_curvature*2
        cavity_length = np.linspace(0.01,max_cavity_length,10000)
        fig_name = str(radius_mirror*2)
        mode_radius = mode_radius_fn(radius_of_curvature,cavity_length,wavelength)
        diffraction_limited_radius_value = diffraction_limited_spot(radius_of_curvature,radius_mirror*2,wavelength)
        mode_radius[diffraction_limited_radius_value > mode_radius] = diffraction_limited_radius_value
        mode_radius_mirror = mode_radius_mirror_fn(cavity_length,mode_radius,wavelength)
        cross_sectional_area_value = cross_sectional_area(mode_radius)
        finesse = finesse_calculations(cavity_length,radius_mirror,mode_radius_mirror)
        cooperativity_values = cooperativity(finesse,wavevector(wavelength),mode_radius)
        # you can remove the quotes and look at the mode waist as a function of cavity length if you would like but i left this commented out
        #plots the mode waist and the waist at the mirror surface
        """
        plt.figure()
        plt.plot(cavity_length, mode_radius,cavity_length,mode_radius_mirror)
        #plt.ylim([0,15])
        plt.show()
        #plots the efficiency of mode collection as a function of mirror diameter
        
        plt.figure()
        plt.plot(cavity_length,efficiency(mode_radius_mirror, core_radius,radius_of_curvature,wavelength),color='black')
        plt.xlim([000,1500])
        plt.title("Coupling Efficiency for Fiber Cavities")
        plt.xlabel(r"Cavity Length ($\mu$m)")
        plt.ylabel("Coupling Efficiency")
        if savefig:
            plt.savefig("efficiency_simulations/"+fig_name+"_mirror_diameter_"+str(radius_of_curvature)+"_ROC_efficiency_curve.png")
            plt.clf()
        else:
            plt.show()
        """
        #plt.figure()
        
        
        plt.plot(cavity_length,finesse,color='black')
        plt.xlim([000,max_cavity_length])
        #plt.hlines(5000, 0,max_cavity_length,color = 'red', linestyle='dashed')
        plt.title("Finesse of Fiber Fabry Perot Cavities")
        plt.xlabel(r"Cavity Length ($\mu$m)")
        plt.ylabel("Finesse")
        if savefig:
            plt.savefig("F_Plots/"+fig_name+"_mirror_diameter_"+str(radius_of_curvature)+"_ROC_Finesse_curve.png")
            plt.clf()
        else:
            plt.show()
        
        plt.figure()
        #finesse = finesse_calculations(cavity_length,radius_mirror,mode_radius_mirror)
        plt.plot(cavity_length,cross_sectional_area_value,color='black')
        plt.xlim([0,max_cavity_length])
        plt.ylim([0,60])
        #plt.hlines(5000, 0,max_cavity_length,color = 'red', linestyle='dashed')
        plt.title("Cross Sectional Area of Fiber Fabry Perot Cavities")
        plt.xlabel(r"Cavity Length ($\mu$m)")
        plt.ylabel(r"Cross Sectional Area ($\mu$m$^2$)")
        if savefig:
            plt.savefig("A_Plots/"+fig_name+"_mirror_diameter_"+str(radius_of_curvature)+"_ROC_Cross_Sectional_Area.png")
            plt.close()
        else:
            plt.show()

        plt.figure()
        plt.plot(cavity_length,finesse/cross_sectional_area_value,color='black')
        plt.xlim([0,max_cavity_length])
        plt.ylim([0,4000])
        #plt.hlines(5000, 0,max_cavity_length,color = 'red', linestyle='dashed')
        plt.title("F/A of Fiber Fabry Perot Cavities")
        plt.xlabel(r"Cavity Length ($\mu$m)")
        plt.ylabel("Finesse/Cross Sectional Area (1/$\mu$m$^2$)")
        if savefig:
            plt.savefig("FA_Plots/"+fig_name+"_mirror_diameter_"+str(radius_of_curvature)+"_ROC_FA_Curve.png")
            plt.close()
        else:
            plt.show()
        plt.figure()
        plt.plot(cavity_length,cooperativity_values,color='black')
        plt.xlim([0,max_cavity_length])
        plt.ylim([0,60])
        #plt.hlines(5000, 0,max_cavity_length,color = 'red', linestyle='dashed')
        plt.title("Cooperativity")
        plt.xlabel(r"Cavity Length ($\mu$m)")
        plt.ylabel("Cooperativity")
        if savefig:
            #plt.savefig("FA_Plots/"+fig_name+"_mirror_diameter_"+str(radius_of_curvature)+"_ROC_FA_Curve.png")
            plt.close()
        else:
            plt.show()
        
