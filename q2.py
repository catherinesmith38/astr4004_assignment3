#Imports
from astroquery.gaia import Gaia


#------------------------------------Task 1------------------------------------
# Find all stars within 1 degree of Ruprecht 147 that are brighter than G mag 14
# cross match with 2MASS catalogue
ra_cen = 289.074 # degrees
dec_cen = -16.323 # degrees

# Remove row limit 
Gaia.ROW_LIMIT = 5000

# Define ADQL query
query = f"""
SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag 
FROM gaiadr3.gaia_source AS g
WHERE DISTANCE({ra_cen}, {dec_cen}, g.ra, g.dec) < 1
AND phot_g_mean_mag < 14
"""

# Cross match query something like this:
#FROM gaiadr3.gaia_source AS g
#JOIN gaiadr3.tmass_psc_xsc_neighbourhood AS tm ON g.source_id = tm.source_id
#JOIN gaiadr1.tmass_original_valid AS tm_orig ON tm.original_ext_source_id = tm_orig.tmass_oid

# Execute the query
job = Gaia.launch_job_async(query)
results = job.get_results()


print(f"Number of stars found: {len(results)}")




#------------------------------------Task 2------------------------------------


