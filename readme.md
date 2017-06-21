# Looking at SN Ia Local Environments with SDSS Data

[![Build Status](https://travis-ci.com/benjaminrose/SNIa-Local-Environments.svg?token=4zcThx8qWuKzVAuCdesD&branch=master)](https://travis-ci.com/benjaminrose/SNIa-Local-Environments) [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 

This project has the goal of calculating the average age of the local environment where a SNIa occurred and seeing if stellar populations age produces a noticeable bias in the Hubble diagram residuals. 

## Code

To reproduce this work, you can follow the `install` section of .travis.yml. I would recommend following FSPS's instructions rather then using the `.travis.Makefile`. The required packages are in the `requirements.txt` file and work nicely with `pip install -r requirements.txt`. 

My research was conducted with python 3.5.2, numpy 1.11.2, scipy 0.18.1, emcee 2.2.1, astropy 1.3, pandas 0.20.1, docpot 0.6.2, FSPS commit [ae31b2f](https://github.com/cconroy20/fsps/commit/ae31b2f63d865354ce944e5c22eba6e93e01e67d), python-fsps commit [6b775a4](https://github.com/dfm/python-fsps/commit/6b775a46cb1cceac145cf08f234f52e04385f001), but should compile with the latest versions of each of these, as tested by Travis-CI.


## Paper

This is building on the work from [this poster](https://ui.adsabs.harvard.edu/#abs/2017AAS...22943402R/abstract) from AAS 229 in January 2017.
