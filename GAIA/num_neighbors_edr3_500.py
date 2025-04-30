'''
This counts the number of neighbors, as defined in Section 2.1, for each star in the search sample.  
Need to run this before making the binary catalog.
MODIFIED VERSION FOR 500-ROW DATASET
'''
from astropy.table import Table
import multiprocessing, psutil
from sklearn.neighbors import BallTree
import numpy as np
from find_binaries_edr3 import duplicates_msk, unique_value_msk, fetch_table_element, get_delta_mu_and_sigma
 
# Use 500-row dataset
tab = Table.read('edr3_parallax_snr5_goodG_500.fits') # 500 sources

size_max_pc = 5 # max projected separation out to which to search
dispersion_max_kms = 5 # max velocity difference in kms
ra, dec, pmra, pmdec, parallax, parallax_error, pmra_error, pmdec_error, G = fetch_table_element(['ra', 'dec', 'pmra', 'pmdec', 'parallax', 'parallax_error', 'pmra_error', 'pmdec_error', 'phot_g_mean_mag'], tab )

s_max_cluster = 206265*size_max_pc
theta_max_radians = s_max_cluster/(1000/parallax)/3600 * np.pi/180
coords = np.vstack([dec*np.pi/180, ra*np.pi/180,]).T
tree = BallTree(coords[G < 18], leaf_size = 10, metric = 'haversine') # build tree of all stars brighter than 18

# data for stars brighter than G = 18
ra_b, dec_b, pmra_b, pmdec_b, parallax_b, parallax_error_b, pmra_error_b, pmdec_error_b, G_b = ra[G < 18], dec[G < 18], pmra[G < 18], pmdec[G < 18], parallax[G < 18], parallax_error[G < 18], pmra_error[G < 18], pmdec_error[G < 18], G[G < 18]

# For 500-row dataset, process in batches
sigma_cut = 2 # how many sigma tolerance

print(f"Processing {len(coords)} stars...")

# Create the array to hold neighbor counts
N_neighbors = np.zeros(len(coords))

# Process 50 stars at a time to show progress
batch_size = 50

for batch_start in range(0, len(coords), batch_size):
    batch_end = min(batch_start + batch_size, len(coords))
    print(f"Processing batch {batch_start//batch_size + 1} of {(len(coords)-1)//batch_size + 1} (stars {batch_start+1}-{batch_end})")
    
    for i in range(batch_start, batch_end):
        # Find companions and their distances
        indices, dists = tree.query_radius([coords[i]], r=theta_max_radians[i], return_distance=True)
        indices, dists = indices[0], dists[0]  # Unpack the nested arrays
        
        # Calculate angular separations
        thetas_arcsec = dists*180/np.pi*3600
        
        # Calculate parallax differences in sigma units
        d_par_over_sigma = np.abs(parallax[i] - parallax_b[indices])/np.sqrt(parallax_error[i]**2 + parallax_error_b[indices]**2)
        
        # Calculate proper motion differences
        delta_mu, sigma_delta_mu = get_delta_mu_and_sigma(
            pmra1=pmra[i], pmdec1=pmdec[i], 
            pmra2=pmra_b[indices], pmdec2=pmdec_b[indices], 
            pmra_error1=pmra_error[i], pmdec_error1=pmdec_error[i], 
            pmra_error2=pmra_error_b[indices], pmdec_error2=pmdec_error_b[indices]
        )
        
        # Maximum allowed proper motion difference based on physical velocity
        mu_max = 0.21095*dispersion_max_kms*parallax[i]
        
        # Count neighbors that satisfy all criteria
        neighbors = (delta_mu < mu_max + sigma_cut*sigma_delta_mu) & (d_par_over_sigma < sigma_cut) & (thetas_arcsec > 1e-3)
        N_neighbors[i] = np.sum(neighbors)

# save these for later
np.savez('neighbor_counts_edr3_500.npz', source_id=fetch_table_element('source_id', tab), N_neighbors=N_neighbors)
print("Saved neighbor counts to neighbor_counts_edr3_500.npz")

# Print summary of neighbor counts
print(f"Stars with >0 neighbors: {np.sum(N_neighbors > 0)}")
print(f"Stars with >30 neighbors (crowded): {np.sum(N_neighbors > 30)}") 
