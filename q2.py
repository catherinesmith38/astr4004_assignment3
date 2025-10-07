#Imports
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits


#------------------------------------Task 1------------------------------------
# Find all stars within 1 degree of Ruprecht 147 that are brighter than G mag 14 and cross match with 2MASS catalogue
ra_cen = 289.074 # degrees
dec_cen = -16.323 # degrees

# Remove row limit 
Gaia.ROW_LIMIT = -1

# Define ADQL query
query = f"""
SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.bp_rp, g.parallax, tm.j_m, tm.h_m, tm.ks_m, tm.ph_qual
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm ON g.source_id = xm.source_id
JOIN gaiadr1.tmass_original_valid AS tm ON xm.clean_tmass_psc_xsc_oid = tm.tmass_oid
WHERE DISTANCE({ra_cen}, {dec_cen}, g.ra, g.dec) < 1
AND phot_g_mean_mag < 14
"""


try:
    # Execute the query
    job = Gaia.launch_job_async(query)
    results = job.get_results()
except:
    # If the query fails use the already downloaded gaia sample
    results = fits.open('data/q2_gaia_sample')[1].data

#------------------------------------Task 2------------------------------------
# Print the number of stars from this query 
print(f"Number of stars found: {len(results)}")


#------------------------------------Task 3------------------------------------
# Identify stars with bad 2MASS photometry
bad_ph = results['ph_qual'] != 'AAA'


#------------------------------------Task 4------------------------------------
# Identify stars with negative parallax
neg_par = results['parallax'] < 0


#------------------------------------Task 5------------------------------------
# Apply quality cuts and print the number of stars remaining
good_stars = results[~bad_ph & ~neg_par]
print(f"Number of good quality stars: {len(good_stars)}")


#------------------------------------Task 6------------------------------------
# Plot CMD from both Gaia and 2MASS

# Calculate absolute G magnitude
distance_pc = 1000 / good_stars['parallax']  # Distance in parsecs
g_abs = good_stars['phot_g_mean_mag'] - 5 * (np.log10(distance_pc) - 1)

# Set up subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Gaia CMD
ax1.scatter(good_stars['bp_rp'], g_abs, s=5, color='blue', alpha=0.5)
ax1.set_xlabel('Bp - Rp (mag)')
ax1.set_ylabel('G (mag)')
ax1.invert_yaxis()
ax1.set_title('Gaia CMD')

# 2MASS CMD
ax2.scatter(good_stars['j_m'] - good_stars['ks_m'], good_stars['h_m'], s=5, color='blue', alpha=0.5)
ax2.set_xlabel('J - Ks (mag)')
ax2.set_ylabel('h (mag)')
ax2.invert_yaxis()
ax2.set_title('2MASS CMD')

fig.suptitle('Color-Magnitude Diagrams for Ruprecht 147')

plt.tight_layout()

#--------------------------------Task 7------------------------------------
# Save the plots

# Make figures directory if one does not already exist
if not os.path.exists('figures'):
    os.makedirs('figures')

# Save figure
plt.savefig('figures/cmds_R147.png', dpi=200)

#--------------------------------Task 8------------------------------------
# Give recommendations for observations 

print('The 2dF fibre positioner on the HERMES spectrograph has 392 science fibres and covers 2 square degrees on the sky')
print('There is a minimum separation of 30 arcseconds between fibres, but the spectrograph can cover most of the defined cone search (3.14 square degrees). ')
print(f'To observed all {len(good_stars)} stars with good photometry, a minimum of {len(good_stars)//392 + 1} pointings are required')