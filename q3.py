#Imports
from astropy.io import fits
import os
from urllib.request import urlretrieve
import numpy as np
import matplotlib.pyplot as plt



#-----------------------------Task 1-----------------------------
# download data and make plot of radial metallictity gradient 

file_path = 'https://github.com/svenbuder/astr4004_2025_week8/blob/main/data/nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits'
file_name = 'nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits'

# if not already downloaded, download the file into data folder
if not os.path.exists('data'):
    # create data directory if it does not exist
    os.makedirs('data')

if not os.path.exists(f'data/{file_name}'):
    # download to data directory 
    urlretrieve(file_path, f'data/{file_name}')

# open the fits file and extract the data
data = fits.open(f'data/{file_name}', ignore_missing_simple=True)[1].data

# extract data 
r_gal = np.sqrt(data['x']**2 + data['y']**2 + data['z']**2)  
a_o = data['a_o']  


# Make plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))


ax1.hist2d(r_gal, a_o, bins=50, cmap='plasma', cmin=1)
ax1.set_xlabel('Galactocentric Radius (kpc)')
ax1.set_ylabel('Gas Phase Metallicity')
#ax1.set_xscale('log')

# linear fit
fit, cov = np.polyfit(r_gal, a_o, 1, cov=True)
x_fit = np.linspace(min(r_gal), max(r_gal), 100)
y_fit = np.polyval(fit, x_fit)
ax1.plot(x_fit, y_fit, color='red', label=f'Fit: y={fit[0]:.2f}x + {fit[1]:.2f}')
ax1.legend(loc='lower left')
ax1.set_title('Radial Metallicity Gradient')

# residuals of fit with log-spaced bins on x-axis
residuals = a_o - np.polyval(fit, r_gal)
ax2.hist2d(r_gal, residuals, bins=50, cmap='plasma', cmin=1)
ax2.axhline(0, color='red', linestyle='--')
ax2.set_xlabel('Galactocentric Radius (kpc)')
ax2.set_ylabel('Residuals')
ax2.set_title('Residuals of Fit')
#ax2.set_xscale('log')
ax2.set_yscale('linear')

plt.tight_layout()

# save to figures 
plt.savefig('figures/rmg_linear_fit.png', dpi=300)


#-----------------------------Task 2-----------------------------
# Quote error on slope and intercept from covariance matrix
slope_err = np.sqrt(cov[0, 0])
intercept_err = np.sqrt(cov[1, 1])

print(f"Slope: {fit[0]:.2e} +/- {slope_err:.2e}")
print(f"Intercept: {fit[1]:.2e} +/- {intercept_err:.2e}")


# -----------------------------Task 3-----------------------------
# Quantify the goodness of fit 

# Root mean square error
rmse = np.sqrt(sum(residuals**2)/len(a_o))
print(f"Root Mean Square Error: {rmse:.2e}")

# R-squared
ss_total = np.sum((a_o - np.mean(a_o))**2)
ss_residual = np.sum(residuals**2)
r_squared = 1 - (ss_residual / ss_total)
print(f"R-squared: {r_squared:.4f}")

