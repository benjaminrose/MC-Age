.. mc-age documentation master file, created by
   sphinx-quickstart on Tue Nov 21 16:48:56 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======================================================
MC-Age --- Estimating Stellar Population Ages with FSPS
=======================================================

About
=====
This project uses FSPS, and python-FSPS to estimate the age of a stellar populations from a photometric SED. The code is available on github_.

.. _github: https://github.com/benjaminrose/SNIa-Local-Environments

Science goal
------------

We check local environment effects on SNIa by looking for correlations between Hubble Residual and the mass weighted average age of the local environment calculated from SDSS Scene Modeling Photometry. We also look at correlations between other variables including SALT2 stretch and color, and host galaxy stellar mass.


.. include:: ../README.rst


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   .. usage
   data-files
   consistency-tests
   script-doc
   API-doc

Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`

Assumptions
------------
* This code assumes you are using SDSS *ugriz* SED's. There are several lines of code that would need to be updated if this assumption is changed. Hopefully they are marked with comments.