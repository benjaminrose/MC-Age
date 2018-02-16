Raw Data
========

This documents where each set of data comes from and how to re-obtain it.

* ``CampellHoltzman.tsv``
    * This is an original data set that combines information from ``SDSS_Photometric_SNe_Ia.fits``
    (Campbell) and ``SMP-photometry`` (Holtzman).
    * It was created via ``Collectin Local SN Data.ipynb``.
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
* ``SDSS_Photometric_SNe_Ia.fits``
    * This is the data released in Campbell 2013. It is the cosmological data
    for the spectroscopic and photometrically classified SN Ia from SDSS.
    * It is originally available at http://www.icg.port.ac.uk/stable/campbelh/SDSS_Photometric_SNe_Ia.fits
