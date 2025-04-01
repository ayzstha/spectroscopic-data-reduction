# Aayush Shrestha | April 16, 2025
# Project for Topic 2: Python Essentials

# Importing necessary libraries
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

# Question 1
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

# Question 2
# (a) Read in the .csv file containing Be star information
filename = "BeStar.csv"

with open(filename, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    be_stars = list(reader)  # Convert to a list of dictionaries

# (b) Print the (HD #, RA, DEC, Magn) values of the 17th-27th Be star entries
print("HD #, RA, DEC, Magn")
for star in be_stars[16:27]:  # Indexing starts from 0, so 17th entry is index 16
    print(star["HD #"], star["RA"], star["DEC"], star["Magn"])

# Question 3
# (a) Read in the FITS file
fits_file = "test_star.fits" 
hdu = fits.open(fits_file)
image_data = hdu[0].data
header = hdu[0].header
hdu.close()

# (b) Plot the image as a grayscale image
plt.figure()
plt.imshow(image_data, cmap='gray', origin='lower', vmin=np.percentile(image_data, 5), vmax=np.percentile(image_data, 95))
plt.colorbar(label='Flux')
plt.title("Star Field Image")

# (c) Use WCS to label axes
wcs = WCS(header)
plt.xlabel("RA (J2000)")
plt.ylabel("Dec (J2000)")
plt.grid(color='white', ls='dotted')
plt.show()

# (d) Read the CSV file and calculate flux means
csv_file = "stars.csv"  # Update with actual file path
stars = pd.read_csv(csv_file)
star_x, star_y = int(stars.iloc[9]['x']), int(stars.iloc[9]['y'])

# Mean flux in +/- 4 pixels around star center
flux_4px = image_data[star_y-4:star_y+5, star_x-4:star_x+5]
mean_flux_4px = np.mean(flux_4px)

# Convert arcseconds to pixels using WCS
arcsec_to_pixel = np.abs(wcs.wcs.cdelt[0]) * 3600
radius_px = int(10 / arcsec_to_pixel)
flux_10arcsec = image_data[star_y-radius_px:star_y+radius_px+1, star_x-radius_px:star_x+radius_px+1]
mean_flux_10arcsec = np.mean(flux_10arcsec)

print(f"Mean flux in +/- 4 pixels: {mean_flux_4px}")
print(f"Mean flux in +/- 10 arcsec: {mean_flux_10arcsec}")

# (e) Create and write new FITS file with log10(flux)
log_flux = np.log10(image_data, where=image_data > 0)  # Avoid log(0)
hdu_new = fits.PrimaryHDU(log_flux, header=header)
hdu_new.writeto("log_star_field.fits", overwrite=True)
print("New FITS file saved as log_star_field.fits")
