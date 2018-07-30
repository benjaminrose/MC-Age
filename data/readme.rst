Detailed Data Set Comments
--------------------------

Detailed comments and documentation on each data set file in the ``data/`` folder.

RA & Dec are in degrees. PetroRad_r is in arcsec. Photometry is in magnitudes.

Main Data
`````````

* ``campbell_local.tsv``
    * An original data set that combines information from ``SDSS_Photometric_SNe_Ia.fits``
    (Campbell) and ``SMP-photometry\`` (Holtzman).
    * Hubble residuals ('hr') use the Malmquist bias corrected distances.
    * It was created via ``Collecting Local SN Data.ipynb`` with the data cuts in the paper: z < 0.5, simga_mag < 1.5 mag, HR < 0.7 mag.
    * N = 104
* ``campbell_global.tsv``
    * The global photometry of the same hosts as ``campbell_local.tsv``
    * Retrieval SEDs in the same way as ``gupta_global.tsv``.
    * N = 104
* ``campbell_formass.csv``
    * TODO
    * a csv version of ``campbell_global.tsv`` 
    * with ``PetroRad_r`` column removed

Calibration Data
````````````````

* ``ellipt_burst_global.tsv``
    * TODO
    * note really SNID's but rather galaxy number. More details are in the ``notes`` header
* ``ellipt_burst_local.tsv``
    * Not needed therefore not there. 
* ``gupta_global.tsv``
    * The objects analyzed in Gupta 2011 but with the same quality cuts I performed on the Campbell/Holtzmann data set.
    * This is the global photometry for the host galaxies of many SDSS-II SN. From DR12 via astroquery 0.3.4, by hand on SkyServer DR13, and Gupta 2011 redshift if none was found before
    * N = 77
* ``Gupta11_table2.???``
    * TODO

``self_consistancey.tsv``
    * no need for RA and Dec in this file. Code skips dust correction for this data set.

H_0 analysis
`````````````

``riess_local.tsv``
``riess_global.tsv``
``riess_formass.csv``

Others
``````

* Table 7
* Table 8
* Table 9

Folders
```````
* ``calibration_sample_lc/``
    * Raw light curves
    * & light curves in the form needed by ``sncosmo`` and will be analyses by  ``get_salt_values.py``
* ``oldruns/``
    * original version of data folder. Used for much of the development. Was removed on 2018-07-27 for this current system.
* ``SDSS - coadded``, ``SDSS - spectra``
    * Not used and are self evident
* ``SMP-photometry``
    * "raw" data for from Holtzman scene modeling.


-------

* ``CampbellHoltzman_formass.csv``
    * This file was used as part of the inputs to ``kcorrect``. 
    * is a subset of ``CampellHoltzman_global.tsv``
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
    * except that these "redshifts" are derived from the median NED cephiad distance and assuming a H_0 = 73.8  km/s/Mpc as fit by Campbell 2013 and used throughout this paper.
* ``SDSS_Photometric_SNe_Ia.fits``
    * This is the data released in Campbell 2013. It is the cosmological data
    for the spectroscopic and photometrically classified SN Ia from SDSS.
    * It is originally available at http://www.icg.port.ac.uk/stable/campbelh/SDSS_Photometric_SNe_Ia.fits
* ``lc/``
    * Raw light curves
    * & light curves in the form needed by ``sncosmo`` and will be analyses by  ``get_salt_values.py``
* ``ellipt_burst_global.tsv``
    * IC 219  host of SN 2005dm   (elliptical
    * NGC 524  host of SN 2000cx  (elliptical)
    * Arp 299 (starburst) host of many CC SN (local is 2010P)
    * M82 (NGC 3034)  (starburst) host of SN 2014J