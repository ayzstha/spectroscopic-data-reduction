# Aayush Shrestha | April 16, 2025

#Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm
import astropy.units as u
from astropy.io import fits
from astropy import wcs
from astropy.io import ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.visualization import make_lupton_rgb
from astropy.modeling import models, fitting
import sys
import csv

# Homework 1
# (a) Create an array X of 100 values from 1 to 100
X = np.arange(1, 101)

# (b) Create array Y that contains X^2
Y = X**2

# (c) Plot X vs. Y
plt.figure()
plt.plot(X, Y, 'r--')  # Red dashed line
plt.xlabel("X Values")
plt.ylabel("Y = X^2")
plt.title("Plot of X vs. Y = X^2")
plt.savefig("plot.pdf")  # Save figure as a .pdf file
plt.show()

# (d) Print Y values that are evenly divisible by 3
print("Y values divisible by 3:")
for y in Y:
    if y % 3 == 0:
        print(y)

# Homework 2
# (a) Read in the .csv file containing Be star information
filename = "BeStar.csv"

with open(filename, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    be_stars = list(reader)  # Convert to a list of dictionaries

# (b) Print the (HD #, RA, DEC, Magn) values of the 17th-27th Be star entries
print("HD #, RA, DEC, Magn")
for star in be_stars[16:27]:  # Indexing starts from 0, so 17th entry is index 16
    print(star["HD #"], star["RA"], star["DEC"], star["Magn"])