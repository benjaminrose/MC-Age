Resources
=========

This folder is the output dumping folder for my projects. For this project these files are too large (over 100 GB) and there for not version controlled (not even with LFS though it was tried at one stage). The final files (used to make the figures in the resulting paper) will be archived somewhere, likely CurateND_. Locally these files are stored in one directory as ``SN{number}_{dataset}_chain.tsv``. The data sets present are:

* circle - fictional objects analyzed for a sanity check
* messier - several messier objects analyzed for a sanity check
* gupta - a reanalysis of Gupta 2011 objects
* campbellG - results for global data of some Campbell 2013 objects
* campbell - results for local data (Holtzman 2008) data of some Campbell 2013 objects

.. _CurateND: https://curate.nd.edu

There is also the resulting stellar mass values stored in ``kcorrect_stellarmass.csv``. These values were generated via kcorrect_ from ``/data/CampbellHoltzman_formass.tsv``.

.. _kcorrect: http://kcorrect.org