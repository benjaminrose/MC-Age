""" main.py -- Looking for correlations between HR and local environments/ages

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-01-19
Python 3
"""

import fsps

"""
# Background

"""

"""
# Goals

"""

"""
# Data Selection

"""

"""
# FSPS Parameters

FSPS has a large number of parameters, many have clear choices, but few need investigation.

## add_igm_absorption

I am going to leave IGM off.

## zcontinuous

We want to use `logzsol` and  let the metallicity of the sample vary a bit (`pmetals`), so we set `zcontinuous = 2`

## cloudy_dust

`True`, dust matters!

## sfh

Looking at Simha 2014, we select the lin-exp model. This get around the systematically over estimate of age of that basic tau-model. Hoping for more accuracy (and justifiably because we are looking at local environments not galaxies as a whole) we will still let t_i vary. (`sfh = 4`)

## compute_vega_mags

The default is `False`, but the examples explicitly set it. Verified below. And passes!

## add_neb_emission
How important is this? I think it is very important for hosts with large positive r-i (> 1). 
"""
import fspsParameters 

# print(fspsParameters.test_compute_vega_mags())
fspsParameters.test_add_neb_emission()