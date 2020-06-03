Stellar Population Age Estimator
================================

.. image:: https://zenodo.org/badge/79574563.svg
   :target: https://zenodo.org/badge/latestdoi/79574563
.. image:: https://travis-ci.org/benjaminrose/MC-Age.svg?branch=master
   :target: https://travis-ci.org/benjaminrose/MC-Age
   :alt: Test status on Travis-CI
.. image:: https://codecov.io/gh/benjaminrose/SNIa-Local-Environments/branch/master/graph/badge.svg?token=sID9V6UFre
	:target: https://codecov.io/gh/benjaminrose/SNIa-Local-Environments
	:alt: Code coverage status via Codecov
.. image:: https://readthedocs.org/projects/mc-age/badge/?version=latest
	:target: https://mc-age.readthedocs.io/en/latest/?badge=latest
	:alt: Documentation Status on Read the Docs
.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
	:target: http://www.astropy.org/
	:alt: This code uses astropy

This project has the goal of calculating the average age of the local environment where a SNIa occurred and seeing if stellar populations age produces a noticeable bias in the Hubble diagram residuals. 

Code
----

Installation
~~~~~~~~~~~~

To reproduce this work, you can follow the ``install`` section of .travis.yml. I would recommend following FSPS's instructions rather then using the ``.travis.Makefile``. The required packages are in the ``requirements.txt`` file and work nicely with ``pip install -r requirements.txt``. To reproduce the plots, you will also need a currently un-aggregated set of packages. Same with the tests and docs.

The dust corrections uses the packages `sfdmap <https://github.com/kbarbary/sfdmap>`_. This requires the download of some fits file dust maps, `hosted on github <https://github.com/kbarbary/sfddata/>`_. This app assumes that these maps are in the folder ``sfddata-master`` at the root of the application. This can be installed by the user, or by running ``bash setup_dust.sh`` in the root directory.

My research was conducted with python 3.5.2, numpy 1.11.2, scipy 0.18.1, emcee 2.2.1, astropy 1.3, pandas 0.20.1, docpot 0.6.2, FSPS commit ae31b2f_, python-fsps commit 6b775a4_, but should "compile" with the latest versions of each of these, as tested by Travis-CI.

.. _ae31b2f: https://github.com/cconroy20/fsps/commit/ae31b2f63d865354ce944e5c22eba6e93e01e67d
.. _6b775a4: https://github.com/dfm/python-fsps/commit/6b775a46cb1cceac145cf08f234f52e04385f001

Usage
~~~~~

A usage message is available at, ::

	python fsps-age.py -h

or as the module doc string in fspsage.py_.

.. _fspsage.py: https://github.com/benjaminrose/SNIa-Local-Environments/blob/master/fspsage.py#L1

Credit
~~~~~~

This code was written by Benjamin Rose. The science was helped by Peter Garnavich. Chris Wotta and Eric Bechter provided feedback on code design issues. Code snippets and external packages are documented in the source code where applied.

Contribution and Licenses
~~~~~~~~~~~~~~~~~~~~~~~~~

I do appreciate scientific verification. If you find an issue or think something should work differently please open an issue or pull request. Specific details and guidelines for contributing will come once interest in this project is shown. If anyone wants to continue this scientific research I would appreciate being contacted about a collaboration.

This software is licensed with the MIT license with the intent to encourage reproducible and open science even though the license is more permissive.

Paper
-----

This is building on the work from this poster_ from AAS 229 in January 2017.

.. _poster: https://ui.adsabs.harvard.edu/#abs/2017AAS...22943402R/abstract

If you use this project or reference the science from it, please cite the forthcoming scientific paper because the code currently does not have a DOI. For now linking to this Github repo is acceptable. DOIs and proper paper citations will come when the project is complete, early Fall 2018.
