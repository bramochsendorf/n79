# n79
These files contain the analysis scripts that support the plots in the 2017 Nature Astronomy letter:

"The star-forming complex LMC-N79 as a future rival to 30 Doradus"

By Bram B. Ochsendorf, Hans Zinnecker, Omnarayani Nayak, John Bally, Margaret Meixner, Olivia C. Jones, Remy Indebetouw, and Mubdi Rahman.

Description of the files:

main.py: 
this script loads in necessary data and calculates star formation and ISM parameters to plot Figure 1, Figure 4, Supplementary Figure 1, and Supplementary Figure 2.

Data used:
- HI (Kim et al., ApJS, 148, 2003)
- CO (Wong et al., ApJS, 197, 16, 2011)
- dust-based molecular gas (Jameson et al., ApJ, 825, 12, 2016)
- HII (Roman-Duval et al., ApJ, 797, 2014)
- Halpha (Gaustadt et al. PASP, 113, 2001)
- 24 micron (Meixner et al., AJ, 132, 2006)

- a list of MYSOs with columns RA, DEC, YSO mass created in build_yso_list.py

aperture_photometry.py: 
this script measures the spectral energy distribution of H72.97-69.39 and plots Figure 3b

build_yso_list.py:
this script takes all the YSO candidates in the LMC, and combines them with the Robitaille et al. 2006 SED fits to extract the stage 1 YSO of M > 8 Msol. It then uses an initial mass function to extimate the total YSO (cluster) mass. The result is used in main.py and plotted in Figure 1c
