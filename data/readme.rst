Raw Data
========

This documents where each set of data comes from and how to re-obtain it.

* ``CampellHoltzman.tsv``
    * This is an original data set that combines information from ``SDSS_Photometric_SNe_Ia.fits``
    (Campbell) and ``SMP-photometry`` (Holtzman).
    * It was created via ``Collecting Local SN Data.ipynb``.
* ``CampbellHoltzman_formass.csv``
    * This file was used as part of the inputs to ``kcorrect``. 
    * is a subset of ``CampellHoltzman_global.tsv``
* ``CampellHoltzman_mb.tsv``
    * This is similar to ``CampellHoltzman.tsv`` but uses Malmquist bias corrected
    distances. No trend with HR is found using this data.
    * It was created via an older copy of ``Collecting Local SN Data.ipynb``.
* ``CampellHoltzman_global.tsv``
	* The global photometry of the same hosts as ``CampellHoltzman.tsv``
	* Retrieval method is the same as ``GobalPhotomtery-Gupta.tsv``.
* ``circlePhotometry.tsv``
	*
* ``galaxyPhotometry.tsv``
	* Am I even using this?
* ``GlobalPhotometry-Gupta.tsv``
	* This is the global photometry for the host galaxies of many SDSS-II SN
From DR12 via astroquery 0.3.4, by hand on SkyServer DR13, and Gupta 2011 redshift if none was found before
* ``Gupta11_table2.txt``
	* The ascii version of table 2 from Gupta 2011. Paper has complete citation.
	* tsv file is the same, I mean to only use the txt file but I appear to be using the tsv once in this whole setup.
* ``Riess2016_calibrators.tsv``
    * These are the galaxies used to calibrate the SN $M_0$ value.
    * They are first reported in Table 1 and used for the middle frame in figure 10.
    * One object is the same as 
    * 101 is M101, 9391 is UGC 9391, and all the rest are NGC objects.
    * The redshifts are from Simbad and the photometry is SDSS modelMag
    * Partial data is below:

'''
101
1365    0.005476    0.00001
1448    0.003895
2442    0.004846
4038    0.005593
5917    0.006472
'''

* ``SDSS_Photometric_SNe_Ia.fits``
    * This is the data released in Campbell 2013. It is the cosmological data
    for the spectroscopic and photometrically classified SN Ia from SDSS.
    * It is originally available at http://www.icg.port.ac.uk/stable/campbelh/SDSS_Photometric_SNe_Ia.fits
