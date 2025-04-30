from astropy.table import Table

# Read the original large file
print("Reading original file...")
tab = Table.read('edr3_parallax_snr5_goodG.fits.gz')

# Extract just the first 500 rows
print(f"Original table has {len(tab)} rows")
small_tab = tab[:500]
print(f"Extracted {len(small_tab)} rows")

# Write to a new file
print("Writing smaller file...")
small_tab.write('edr3_parallax_snr5_goodG_500.fits', format='fits', overwrite=True)
print("Done! Created edr3_parallax_snr5_goodG_500.fits") 
