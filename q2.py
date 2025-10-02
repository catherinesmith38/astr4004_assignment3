#Imports
from astroquery.gaia import Gaia


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
WHERE DISTANCE(289.074, -16.323, g.ra, g.dec) < 1
AND phot_g_mean_mag < 14
"""

#TODO get the cross match working
# Cross match query looks something like this
#FROM gaiadr3.gaia_source AS g
#JOIN gaiadr3.tmass_psc_xsc_neighbourhood AS tm ON g.source_id = tm.source_id
#JOIN gaiadr1.tmass_original_valid AS tm_orig ON tm.original_ext_source_id = tm_orig.tmass_oid

# Execute the query
job = Gaia.launch_job_async(query)
results = job.get_results()


#------------------------------------Task 2------------------------------------


print(len(results))

