Raw Data
========

This documents each data set.

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
* ``final_mags.txt``
    * Output from Peter's magnitude calculations.
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
    * 101 is M101, 9391 is UGC 9391, and all the rest are NGC objects.
    * The redshifts are from Simbad
    * The photometry was done by Peter emailed to me on 2018-06-27 with the subject "final mags for cepheid hosts"
    * "The uncertainties should be dominated by systematics and not statistical noise given the brightness of these objects. I would say that the relative uncertainty from band to band is 0.03 mag. The over-all uncertainty is larger as a slightly bigger aperture might include more light. But use 0.03 mag and we can inflate the errors on the mass." -- Peter
    * Some objects are outside of SDSS footprint, but sadly PanSTARS did not observe in the *u* band:

'''
1365    0.005476    0.00001
1448    0.003895
2442    0.004846
4038    0.005593
5917    0.006472
'''
* ``Riess2016_calibrators_local.tsv``
    * Same as ``Riess2016_calibrators.tsv`` but for 1.5kpc radii apertures around the location of their SN.
* ``Riess2016_formass.csv``
    * same as ``Riess2016_calibrators.tsv`` but saved as a csv rather than tsv.
* ``SDSS_Photometric_SNe_Ia.fits``
    * This is the data released in Campbell 2013. It is the cosmological data
    for the spectroscopic and photometrically classified SN Ia from SDSS.
    * It is originally available at http://www.icg.port.ac.uk/stable/campbelh/SDSS_Photometric_SNe_Ia.fits
