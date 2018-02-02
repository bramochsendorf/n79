# n79
These files contain the analysis scripts that support the plots in the 2017 Nature Astronomy letter:

"The star-forming complex LMC-N79 as a future rival to 30 Doradus"

By Bram B. Ochsendorf, Hans Zinnecker, Omnarayani Nayak, John Bally, Margaret Meixner, Olivia C. Jones, Remy Indebetouw, and Mubdi Rahman.

Data used:
- HI (Kim et al., ApJS, 148, 2003)
- CO (Wong et al., ApJS, 197, 16, 2011)
- 'dust-based' molecular gas (Jameson et al., ApJ, 825, 12, 2016)
- HII (Roman-Duval et al., ApJ, 797, 2014)
- Halpha (Gaustadt et al. PASP, 113, 2001)
- 24 micron (Meixner et al., AJ, 132, 2006)
- a list of MYSOs with columns RA, DEC, YSO mass created in build_yso_list.py

Description of the main analysis scripts:

main.py: 
this script loads in data and calculates star formation and ISM parameters to plot Figure 1, Figure 4, Supplementary Figure 1, and Supplementary Figure 2.

plots.py:
script that is called in main.py and makes the aforementioned plots.

aperture_photometry.py: 
this script measures the spectral energy distribution of H72.97-69.39 and plots Figure 3b.

build_yso_list.py:
This script takes all the YSO candidates in the LMC (output/yso_list/catalogue_merged_yso_clumps_uplimit.txt), and combines them with the Robitaille et al. 2006 SED fits (output/yso_list/merged_ysos_clumps_nopacs100_uplimit.txt) to extract the Stage 1 YSOs of M > 8 Msol. It then uses an initial mass function to extimate the total YSO (cluster) mass of each source. The result (sfr_radec_myso.txt) is used in main.py and used in Figure 1c.
