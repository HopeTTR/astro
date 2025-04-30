'''
This makes the binary catalog that accompanies El-Badry et al. 2021.  
MODIFIED VERSION FOR 500-ROW DATASET
'''

from astropy.table import Table
import multiprocessing, psutil
from sklearn.neighbors import BallTree
import numpy as np
from math import log10
import os

# Define helper functions first - these are needed by num_neighbors_edr3.py

def fetch_table_element(colname, table):
    '''
    avoid table['col'].data vs table['col'].data.data problems with masked arrays in astropy tables 
    '''
    if type(colname) == str:
        if colname in table.colnames:
            if type(table[colname].data.data) == memoryview:
                dat_ = table[colname].data
            else:
                dat_ = table[colname].data.data
        else:
            # If column doesn't exist, return None
            dat_ = None
    elif type(colname) == list:
        dat_ = []
        for col in colname:
            dat_.append(fetch_table_element(col, table))
    return dat_

def find_y_in_x(x, y):
    '''
    x and y are arrays. find the indices in x where the array is equal to values in y. 
    '''
    index = np.argsort(x)
    sorted_x = x[index]
    sorted_index = np.searchsorted(sorted_x, y)
    yindex = np.take(index, sorted_index, mode = "clip")
    mask = x[yindex] == y
    return yindex[mask]

def duplicates_msk(Array):
    '''
    Finds duplicate values of an array.
    Uses masking. Here Array should be a numpy array of ints. Behavior is undefined if array is of 
    lists, sets, etc. 
    '''
    Array = np.array(Array)
    m = np.zeros_like(Array, dtype=bool)
    m[np.unique(Array, return_index=True)[1]] = True
    return ~m
         
def unique_value_msk(x):
    '''
    returns true only if that value appears only once in the array
    '''
    uniq, counts = np.unique(x, return_counts=True)  
    unique_vals = uniq[counts == 1]
    w = find_y_in_x(x, unique_vals)
    is_unique = np.zeros(len(x), dtype = bool)
    is_unique[w] = True
    return is_unique
                    
def get_distance_arcsec(ra1, dec1, ra2, dec2):
    '''
    angular separations. coords are assumed to be in degrees
    '''
    ra_rad1, dec_rad1 = ra1*np.pi/180, dec1 * np.pi/180
    ra_rad2, dec_rad2 = ra2*np.pi/180, dec2 * np.pi/180
    d_ra, d_dec = ra_rad1 - ra_rad2, dec_rad1 - dec_rad2
    
    d_theta = 2*np.arcsin(np.sqrt(np.sin(0.5*d_dec)**2 + np.cos(dec_rad1)*np.cos(dec_rad2)*np.sin(0.5*d_ra)**2))
    d_theta_deg = 180/np.pi*d_theta
    d_theta_arcsec = d_theta_deg * 3600
    return d_theta_arcsec
     
def get_delta_mu_and_sigma(pmra1, pmdec1, pmra2, pmdec2, pmra_error1, 
    pmdec_error1, pmra_error2, pmdec_error2):
    '''
    Uses standard uncertainty propagation 
    Equations 4-5 of the paper. 

    assume that "1" is a float and "2" is an array
    '''
    delt_alpha, delt_delta = (pmra1 - pmra2)**2, (pmdec1 - pmdec2)**2
    delta_mu2 = delt_alpha + delt_delta
    
    try:
        lenn = len(pmra2) # checks whether pmra2 is an array 
        m = delta_mu2 == 0
        sigma2_delta_mu = np.zeros(len(pmra2))
        if np.sum(m):
            sigma2_delta_mu[m] = (pmra_error1**2 + pmra_error2[m]**2) + (pmdec_error1**2 + pmdec_error2[m]**2)
        if np.sum(~m):
            sigma2_delta_mu[~m] = ((pmra_error1**2 + pmra_error2[~m]**2) * (delt_alpha[~m]) + \
                (pmdec_error1**2 + pmdec_error2[~m]**2)*delt_delta[~m])/delta_mu2[~m]
    except: # pmra2 is a float 
        if delta_mu2 == 0:
            sigma2_delta_mu = pmra_error1**2 + pmra_error2**2 + pmdec_error1**2 + pmdec_error2**2
        else:
            sigma2_delta_mu = ((pmra_error1**2 + pmra_error2**2) * (delt_alpha) + \
                (pmdec_error1**2 + pmdec_error2**2)*delt_delta)/delta_mu2

    return np.sqrt(delta_mu2), np.sqrt(sigma2_delta_mu)

# Only run the main binary catalog creation if this script is executed directly
if __name__ == "__main__":
    # Constants for the binary search
    parallax_sigma_limit = 3 # only accept pair with parallaxes within 3 sigma of each other
    theta_arcsec_min = 4 # limit below which we'll accept parallaxes within 6 sigma of each other. 
    
    # Use the 500-row dataset
    tab = Table.read('edr3_parallax_snr5_goodG_500.fits')  # 500 elements

    # remove stars that have too many neighbors, as defined in section 2.1 
    # if you want to look at the "initial candidates" sample, including clusters, comment this out. 
    if os.path.exists('neighbor_counts_edr3_500.npz'):
        tmp = np.load('neighbor_counts_edr3_500.npz')
        crowded = np.in1d(fetch_table_element('source_id', tab), tmp['source_id'][tmp['N_neighbors'] > 30]); tmp.close()
        tab = tab[~crowded] # remove crowded stars
        print(f"After removing crowded stars, {len(tab)} stars remain")


    source_id, ra, dec, G, parallax, parallax_error, pmra, pmdec, pmra_error, pmdec_error  = fetch_table_element(['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error'], tab)

    s_max_au = 3600*180/np.pi # 206265 au = 1 pc
    theta_max_radians = s_max_au/(1000/parallax)/3600 * np.pi/180 # angular separation corresponding to s = 1 pc
    coords = np.vstack([ dec*np.pi/180, ra*np.pi/180]).T
    tree = BallTree(coords, leaf_size = 20, metric = 'haversine')
    print('built tree') 

    # Process stars in batches (for 500 row dataset)
    print("Searching for binary candidates...")
    star1s, star2s = [], []
    
    # Process in batches for better progress tracking
    batch_size = 50
    
    for batch_start in range(0, len(coords), batch_size):
        batch_end = min(batch_start + batch_size, len(coords))
        print(f"Processing batch {batch_start//batch_size + 1} of {(len(coords)-1)//batch_size + 1} (stars {batch_start+1}-{batch_end})")
        
        batch_companions = 0
        
        for i in range(batch_start, batch_end):
            # Get neighbors within separation limit
            idxs, dists = tree.query_radius([coords[i]], r = theta_max_radians[i], return_distance = True)
            idxs, dists = idxs[0], dists[0]  # Unpack the nested arrays
            
            # Calculate angular separations
            thetas_arcsec = dists*180/np.pi*3600
            
            # Get brighter parallax
            brighter_parallax = np.copy(parallax[idxs]) 
            brighter_parallax[G[idxs] > G[i]] = parallax[i]  # parallax of the brighter component
        
            # Calculate parallax difference in sigma units
            d_par_over_sigma = np.abs(parallax[i] - parallax[idxs])/np.sqrt(parallax_error[i]**2 + parallax_error[idxs]**2)
            
            # Calculate proper motion differences
            delta_mu, sigma_delta_mu = get_delta_mu_and_sigma(
                pmra1 = pmra[i], pmdec1 = pmdec[i], 
                pmra2 = pmra[idxs], pmdec2 = pmdec[idxs], 
                pmra_error1 = pmra_error[i], pmdec_error1 = pmdec_error[i], 
                pmra_error2 = pmra_error[idxs], pmdec_error2 = pmdec_error[idxs]
            )
            
            # Avoid divided-by-zero warnings when calculating delta_mu_orbit for theta = 0
            delta_mu_orbit = np.zeros(len(thetas_arcsec))
            mm = thetas_arcsec == 0
            delta_mu_orbit[mm] = 1e9  # effectively infinite limit for self-pairing
            delta_mu_orbit[~mm] = 0.44428*brighter_parallax[~mm]**(3/2)*thetas_arcsec[~mm]**(-1/2)
            
            # Calculate physical separation in AU
            sep_AU = 1000/brighter_parallax * thetas_arcsec
            
            # Parallax difference limit varies with angular separation
            max_parallax_diff = np.ones(len(thetas_arcsec))*parallax_sigma_limit
            max_parallax_diff[thetas_arcsec < theta_arcsec_min] = 2*parallax_sigma_limit
            
            # Apply all criteria to find binary candidates
            m = (d_par_over_sigma < max_parallax_diff) & \
                (delta_mu < delta_mu_orbit + 2*sigma_delta_mu) & \
                (thetas_arcsec > 0.001) & \
                (sep_AU < s_max_au)
                
            if np.sum(m):
                batch_companions += np.sum(m)
                for k in range(np.sum(m)):
                    star1s.append(i)
                    star2s.append(idxs[m][k])
        
        print(f"Batch complete. Found {batch_companions} possible companions in this batch.")
    
    star1s = np.array(star1s)
    star2s = np.array(star2s)
    print(f'Total length of preliminary catalog is {len(star1s)}')

    # make a new table. each row corresponds to a different pair.
    new_cat = Table()    
    for col in tab.colnames:
        new_cat[col+'1'] = tab[col][star1s]
        new_cat[col+'2'] = tab[col][star2s]

    # remove duplicates (pairs where star 1 and star 2 are switched)
    sid1, sid2 = fetch_table_element(['source_id1', 'source_id2'], new_cat)
    joint_ids = np.vstack([sid1, sid2]).T
    sorted_joint = np.sort(joint_ids, axis=1)
    joint_1d = np.core.defchararray.add(sorted_joint.T[0].astype(str), sorted_joint.T[1].astype(str))
    dups = duplicates_msk(joint_1d)   
    new_cat = new_cat[~dups]
    print(f'After finding {np.sum(dups)} exact duplicates, the new length is {len(new_cat)}')

    # make the brighter star the star "1" and the fainter star "2"
    G1, G2 = fetch_table_element(['phot_g_mean_mag1', 'phot_g_mean_mag2'], new_cat)
    switch = G1 > G2
    colnames =  [ c[:-1] for c in new_cat.colnames if c[-1]=='1'] 
    for col in colnames:
        new_cat[col+'1'][switch], new_cat[col+'2'][switch] = new_cat[col+'2'][switch], new_cat[col+'1'][switch]

    # calculate angular and physical separations. 
    ra1, dec1, ra2, dec2, parallax1, id1, id2 = fetch_table_element(['ra1', 'dec1', 'ra2', 'dec2', 'parallax1', 'source_id1', 'source_id2'], new_cat)
    theta_arcsec = get_distance_arcsec(ra1 = ra1, dec1 = dec1, ra2 = ra2, dec2 = dec2)
    new_cat['pairdistance'] = theta_arcsec/3600
    new_cat['sep_AU'] = 1000/parallax1 * theta_arcsec

    # remove triples. first get rid of cases where id1 or id2 have duplicates; i.e. one star with two companions 
    new_cat = new_cat[unique_value_msk(id1) & unique_value_msk(id2)]

    # now get cases where elements of id1 are in id2 or vice versa.
    id1, id2 = fetch_table_element(['source_id1', 'source_id2'], new_cat)
    new_cat = new_cat[~ (np.in1d(id1, id2) | np.in1d(id2, id1))]
    print(f'After removing triples, there are {len(new_cat)} pairs')

    # Check if bp_rp columns exist in the dataset
    has_color = 'bp_rp' in tab.colnames
    
    if has_color:
        # Classify star types if color information is available
        print("Color information found. Classifying binary types...")
        bp_rp1, G1, bp_rp2, G2, parallax = fetch_table_element(['bp_rp1', 'phot_g_mean_mag1', 'bp_rp2', 'phot_g_mean_mag2', 'parallax1'], new_cat)
        Mg1 = G1 + 5 * log10(parallax / 100) 
        Mg2 = G2 + 5 * log10(parallax / 100) 

        col1, col2 = np.isfinite(bp_rp1), np.isfinite(bp_rp2)  # whether or not each component has a color

        wd1 = (Mg1 > (3.25*bp_rp1 + 9.625)) & col1
        wd2 = (Mg2 > (3.25*bp_rp2 + 9.625)) & col2
        wdms = (wd1 & ~wd2) | (~wd1 & wd2) & col1 & col2
        wdwd = wd1 & wd2  & col1 & col2
        msms = ~wd1 & ~wd2 & col1 & col2
        wd__ = wd1 & col1 & ~col2
        __wd = wd2 & col2 & ~col1
        ms__ = (~wd1 & col1) & ~col2
        __ms = (~wd2 & col2) & ~col1
        ____ = ~col1 & ~col2

        binary_type = np.array(len(new_cat)*['XXXX'])
        binary_type[wdwd] = 'WDWD'
        binary_type[wdms] = 'WDMS'
        binary_type[msms] = 'MSMS'
        binary_type[wd__ | __wd] = 'WD??'
        binary_type[ms__ | __ms] = 'MS??'
        binary_type[____] = '????'
    else:
        # Skip the star type classification since we're missing the bp_rp columns
        print("Note: Star type classification skipped as bp_rp columns are missing in the dataset")
        binary_type = np.array(len(new_cat)*['UNKNOWN'])
    
    new_cat['binary_type'] = binary_type

    # write out the catalog in both FITS and CSV formats
    new_cat.write('binary_catalog_500.fits', format='fits', overwrite=True)
    new_cat.write('binary_catalog_500.csv', format='csv', overwrite=True)
    
    print("Created binary catalog: binary_catalog_500.fits and binary_catalog_500.csv")
    
    # Print a sample and summary of the results
    print("\nSample of the binary catalog (first few rows):")
    for i in range(min(5, len(new_cat))):
        print(f"Pair {i+1}: Separation = {new_cat['sep_AU'][i]:.1f} AU, Angular separation = {new_cat['pairdistance'][i]*3600:.2f} arcsec")
    
    print(f"\nSummary: Found {len(new_cat)} binary candidates from {len(tab)} stars") 
