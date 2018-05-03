
# coding: utf-8

# In[5]:

import numpy as np
import astropy.units as u

def Read(filename):
    # Reads the data from a file
    file = open(filename, 'r') # Opens file
    line1 = file.readline()
    label, value = line1.split()
    time = float(value) * u.Myr ## This is the timestep the file represents
    line2 = file.readline()
    label2, value2 = line2.split()
    total = float(value2) # This is the total number of data points in the file
    file.close()
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3) # The nonheader data
    return time, total, data


# In[ ]:



