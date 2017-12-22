Consistency Tests
=================

Circle Test
-----------

The first test is to make sure the MCMC can get out the same result that ``FSPS`` was given. Seeing the flow of the values 

+-------------------------+----------------------+-------------------+----------------------+
|Input 1 (to FSPS)        | Output 1 (from FSPS) | Input 2 (to spae) | Output 2 (from spea) |
+=========================+======================+===================+======================+
|SFH parameters, redshift | SED                  | SED, redshift     | SFH parameters       |
+-------------------------+----------------------+-------------------+----------------------+

So the output of ``spea`` should be the same as the input to ``FSPS``.


Setting and generating the values in ``data/circlePhotometry.tsv``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lets generate the needed stellar population. This should be the same as the default in calcualteAge.pyL427_.

.. _calcualteAge.pyL427: https://github.com/benjaminrose/SNIa-Local-Environments/blob/master/calculateAge.py#L427

.. code:: python

    import fsps
    sp = fsps.StellarPopulation(zcontinuous=2, cloudy_dust=True,
                                add_neb_emission = True, sfh=5)
    sdss_bands = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']

Here are the different input 1, star formation histories in a table

+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|name | logzsol | dust2 | tau | tStart | sfTrans | sfSlope | description                |
+=====+=========+=======+=====+========+=========+=========+============================+
|c1   | -0.5    | 0.1   | 0.5 | 1.5    | 9.0     | -1.0    | old                        |
+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|c2   | -0.5    | 0.1   | 0.5 | 1.5    | 9.0     | 15.0    | sharp burst & young        |
+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|c3   | -0.5    | 0.1   | 7.0 | 3.0    | 10.0    | 15.0    | flat burst & young         |
+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|c4   | -0.5    | 0.1   | 7.0 | 3.0    | 13.0    | 0.0     | flat burst                 |
+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|c5   | -1.5    | 0.1   | 0.5 | 1.5    | 9.0     | -1.0    | metal poor, like c1        |
+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|c6   | -0.5    | 0.8   | 7.0 | 3.0    | 10.0    | 15.0    | dusty, like c2             |
+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|c7   | -0.5    | 0.1   | 0.5 | 1.5    | 6.0     | 15.0    | mostly late linear, like c2|
+-----+---------+-------+-----+--------+---------+---------+----------------------------+
|c8   | -0.5    | 0.1   | 0.1 | 8.0    | 12.0    | 20.0    | very young                 |
+-----+---------+-------+-----+--------+---------+---------+----------------------------+

and as copyable code

.. code:: python

    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.1, 0.5, 1.5, 9.0, -1.0
    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.1, 0.5, 1.5, 9.0, 15.0
    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.1, 7.0, 3.0, 10.0, 15.0
    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.1, 7.0, 3.0, 13.0, 0.0
    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -1.5, 0.1, 0.5, 1.5, 9.0, -1.0
    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.8, 7.0, 3.0, 10.0, 15.0
    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.1, 0.5, 1.5, 6.0, 15.0
    logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.1, 0.1, 8.0, 12.0, 20.0

.. logzsol, dust2, tau, tStart, sfTrans, sfSlope = 0.0, 0.1, 10.0, 6.0, 8.0, 19.0
.. logzsol, dust2, tau, tStart, sfTrans, sfSlope = 0.0, 0.1, 10.0, 9.0, 10.0, 19.0

.. # 2017-08-24 results of Circle 1
.. logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.69, 0.21, 2.63, 3.99, 4.04, 6.82
.. c = -25.70

.. code:: python

    sp.params['logzsol'] = logzsol
    dust1 = 2.0*dust2
    sp.params['dust1'] = dust1
    sp.params['dust2'] = dust2
    sp.params['tau'] = tau
    sp.params['sf_start'] = tStart
    sp.params['sf_trunc'] = sfTrans
    sp.params['sf_slope'] = sfSlope
    sp.get_mags(tage=13.185, redshift=0.05, bands=sdss_bands)

And the results can be found in ``data/circlePhotometry.tsv`` available on GitHub_ with a boost of c = 25

.. _GitHub: https://github.com/benjaminrose/SNIa-Local-Environments/blob/master/data/circlePhotometry.tsv

comment on photometric uncertainties.

SED Check
~~~~~~~~~~

Messier Objects
---------------

The paper coves this well.

Gupta Objects
-------------

The paper covers this well.


.. # these should have gotten: logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.5, 0.1, 0.5, 1.5, 9.0, -1.0 with c = -25
.. # mangitudes = 
.. # logzsol can be as low as -2.5
.. logzsol, dust2, tau, tStart, sfTrans, sfSlope = -2.5, 0.01, 7.17, 7.94, 10.40, -5.24
.. c = -23.48
.. instruments
.. array([ 43.31499986,  42.05981015,  41.76335814,  41.66908196,  41.62473476])
.. # logzsol can be as low as -1.0
.. logzsol, dust2, tau, tStart, sfTrans, sfSlope = -1.0, 0.25, 5.67, 1.94, 4.93, 1.64
.. c = -22.85
.. array([ 42.27927002,  41.4316055 ,  41.23214939,  41.01247384,  40.99344229])
.. # logzsol is fixed at -0.5
.. logzsol, dust2, tau, tStart, sfTrans, sfSlope = -0.51, 0.32, 8.17, 8.42, 10.76, 4.72
.. c = -22.17
.. array([ 41.5268658 ,  40.70263292,  40.54702704,  40.32837067,  40.30277849])