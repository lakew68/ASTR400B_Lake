
# coding: utf-8

# In[1]:

# Project Piece 1
# Velocity dispersion
# William Lake

# import modules                                                                                                       
import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterOfMass2 import CenterOfMass
import math
import matplotlib.pyplot as plt
plt.ioff()

from astropy.constants import G
G = G.to(u.kpc * u.km ** 2 / u.s ** 2 / u.Msun) # Converts the gravitational constant to our desired units


# In[44]:

class VelocityDispersion:
    '''This class aids in finding the velocity dispersions of bulge particles over radius. It will eventually
    do the same over time, held at constant radius'''
    
    def __init__(self, galaxy, snap):
        # This function/class takes the galaxy name and snapshot number as inputs
        
        # The following three lines create the part of the filename that describe the Snap number, and then generate the filename
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]
        self.filename = "HighRes/%s_"%(galaxy) + ilbl + '.txt'
        time, total, data = Read(self.filename) # Read the data
        self.x = np.array(data['x']) * u.kpc # Import the x, y, and z coordinates with correct units
        self.y = np.array(data['y']) * u.kpc
        self.z = np.array(data['z']) * u.kpc
        self.m = np.array(data['m']) # Import the mass data
        self.vx = np.array(data['vx']) * u.km / u.s
        self.vy = np.array(data['vy']) * u.km / u.s
        self.vz = np.array(data['vz']) * u.km / u.s
        self.data = data # Wasn't asked for but useful later
        self.gname = galaxy # Stores galaxy name
        
    def distance(self, x1, x2, y1, y2, z1, z2):
        # This function describes the absolute distance between points (x1,y1,z1) and (x2, y2, z2)
        return np.sqrt(((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2))
    
    def dispersion(self, vArray, VX_COM, VY_COM, VZ_COM):
        # This function finds the velocity dispersion of an array vArray of 3-D velocities (contained within a formatted data set)
        bulgeIndices = np.where(vArray['type'] == 3)
        
        corr_vArray = []
        bulgeVArray = vArray[bulgeIndices]
        
        for i in range(np.shape(bulgeVArray)[0]):
            corr_vArray.append(self.distance(bulgeVArray['vx'][i] * u.km / u.s, VX_COM, bulgeVArray['vy'][i] * u.km / u.s, VY_COM, bulgeVArray['vz'][i] * u.km / u.s, VZ_COM).value)
            
        disp = np.std(corr_vArray)
        return disp
        
    def CumDispersion(self, rArray):
        # Takes as input an array of radii
        # Returns the cumulative velocity dispersion within the stars enclosed
        COM_Object = CenterOfMass(self.filename, 2)
        VX_COM, VY_COM, VZ_COM = COM_Object.COM_V(1.0)
        X_COM, Y_COM, Z_COM = COM_Object.COM_P(1.0) # Finds the center of mass coordinates of the galaxy
        
        typeIndex = np.where(self.data['type'] == 3) # Selects bulge particles
        
        radii = self.distance(self.x, X_COM, self.y, Y_COM, self.z, Z_COM) 
        # Creates an array to hold distances of particles from the COM
        dispersions = np.zeros(len(rArray)) # Creates our result array
        
        for index in range(len(rArray)):
            # Loops over rArray to find enclosed dispersion
            indices = np.where(radii < rArray[index])
            dispersions[index] = self.dispersion(self.data[indices], VX_COM, VY_COM, VZ_COM)
            print(rArray[index])
        
        return dispersions


# In[45]:

MWDisp = VelocityDispersion('MW', 800)

radiusArray = []
i = 30
while i > .4:
    radiusArray.append(i)
    i /= 1.1
radiusArray = np.array(radiusArray) * u.kpc

dispersion = MWDisp.CumDispersion(radiusArray)
print(dispersion)
fig = plt.figure()

plt.plot(radiusArray, dispersion,'b', label = 'MW Cumulative Dispersion')

plt.savefig('MW_800_Dispersion_Plot.png')

plt.close(fig)