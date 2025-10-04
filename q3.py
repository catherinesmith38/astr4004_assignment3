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
slope = fit[0]
intercept = fit[1]

slope_err = np.sqrt(cov[0, 0])
intercept_err = np.sqrt(cov[1, 1])

print(f"Slope: {slope:.2e} +/- {slope_err:.2e}")
print(f"Intercept: {intercept:.2e} +/- {intercept_err:.2e}")


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

# -----------------------------Task 4-----------------------------
# plot 3 planel figure of 2d histograms in the x y plane of median a_o, linear fit and residuals

# collapse a_o measurements in the z direction and take the median value
def median_xy(x, y, variable, bins):
    """
    Parameters
    x: x coordinates of the data to be collapsed.
    y: y coordinates of the data to be collapsed.
    variable: The data to be collapsed.
    bins: number of bins to calculate the median over.

    Output:
    median: matrix containing median values of varible. Size bins x bins.
    (x_bin_edges, y_bin_edges): bin edges in x and y direction for the matrix.
    """

    # define x and y bins 
    x_bin_edges = np.linspace(min(x), max(x), bins+1)
    y_bin_edges = np.linspace(min(y), max(y), bins+1)


    # Collect a_o measurements in these x and y bins
    # initiate matrix 
    median = np.empty(shape=(bins, bins))

    # loop over bins and compute median a_o for all gas particles in each (x, y) cell
    for i in range(bins):
        for j in range(bins):
            # mask of points within the (x,y) bin
            mask = (
                (x >= x_bin_edges[i]) & (x < x_bin_edges[i+1]) &
                (y >= y_bin_edges[j]) & (y < y_bin_edges[j+1])
            )
            # compute median a_o if there are data points in the bin
            if np.any(mask):
                median[i, j] = np.median(variable[mask])
            else:
                # if there is no data set the value to nan
                median[i, j] = np.nan
    return median, (x_bin_edges, y_bin_edges)

# take collapsed median of simulated data
a_o_median, (x_bin_edges, y_bin_edges) = median_xy(data['x'], data['y'], data['a_o'], bins=30)


# Calculate fit 
a_o_fit = slope * np.sqrt(data['x']**2 + data['y']**2 + data['z']**2) + intercept
# Take the collapsed median of fitted data 
a_o_fit_median, (x_bin_edges_fit, y_bin_edges_fit) = median_xy(data['x'], data['y'], a_o_fit, bins=30)

# Find residuals
residuals = data['a_o'] - a_o_fit
# take the collapsed median
residuals_median, (x_res, y_res) = median_xy(data['x'], data['y'], residuals, bins=30)

# set up figure
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(15, 6))

# median simulated a_o
im1 = ax1.imshow(
    a_o_median.T,
    origin='lower',
    extent=[x_bin_edges[0], x_bin_edges[-1], y_bin_edges[0], y_bin_edges[-1]],
    cmap='plasma',
    aspect='equal', 
    vmax=np.nanmax([np.nanmax(a_o_median), np.nanmax(a_o_fit_median)]), 
    vmin=np.nanmin([np.nanmin(a_o_median), np.nanmin(a_o_fit_median)]) # ensure fitted and simulated data are on the same colour scale
)
ax1.set_xlabel('X (kpc)')
ax1.set_ylabel('Y (kpc)')
ax1.set_title('Median Simulated A(O)')
plt.colorbar(im1, ax=ax1, label='A(O)', location='bottom')

# median fit a_o
im2 = ax2.imshow(
    a_o_fit_median.T,
    origin='lower',
    extent=[x_bin_edges_fit[0], x_bin_edges_fit[-1], y_bin_edges_fit[0], y_bin_edges_fit[-1]],
    cmap='plasma',
    aspect='equal', 
    vmax=np.nanmax([np.nanmax(a_o_median), np.nanmax(a_o_fit_median)]), 
    vmin=np.nanmin([np.nanmin(a_o_median), np.nanmin(a_o_fit_median)])
)
ax2.set_xlabel('X (kpc)')
ax2.set_ylabel('Y (kpc)')
ax2.set_title('Median Fit A(O)')
plt.colorbar(im2, ax=ax2, label='Fit A(O)', location='bottom')

# median residuals
im3 = ax3.imshow(
    residuals_median.T,
    origin='lower',
    extent=[x_res[0], x_res[-1], y_res[0], y_res[-1]],
    cmap='PRGn',  # diverging colour map
    aspect='equal',
    vmin=-np.nanmax(np.abs(residuals_median)),
    vmax=np.nanmax(np.abs(residuals_median))   # ensure colour map is centred on zero
)
ax3.set_xlabel('X (kpc)')
ax3.set_ylabel('Y (kpc)')
ax3.set_title('Median Residuals')
plt.colorbar(im3, ax=ax3, label=r'$\Delta$A(O)', location='bottom')

plt.tight_layout()
# save to figures
plt.savefig('figures/histograms_a_o_fit_residuals.png', dpi=300)

