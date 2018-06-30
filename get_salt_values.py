""" get_salt_values.py - gets the SALT parameters for H_0 calibration SN

This script assumes:
* light curve files are at ``data/calibration_sample_lc/'``
* light curve files will be named like ``2011fe.txt``
* the light curve files are formated according to ``sncosmo``
* the light curve files have meta data for ``ra``, ``dec`` and ``z``
* RA and Dec meta data are in degrees in the ICRS coordinate system (e.g., "J2000")
* ``SFD_DIR`` environment variable is set according to https://github.com/kbarbary/sfdmap
* Does not robustly handle paths/directories

* Benjamin Rose
* University of Notre Dame
* benjamin.rose@me.com
* Python 3.6
* 2018-05-30
"""

import sncosmo
import sfdmap
import matplotlib.pyplot as plt

dustmap = sfdmap.SFDMap()    # using ``SFD_DIR`` environment variable
dust = sncosmo.CCM89Dust()

sn_list = ['2011fe', '2009ig', '2002fk', '1995al', '1994ae', '2012ht', '2011by', '1998aq',
           '2012cg', '1981B', '1990N', '2007af', '2013dy', '2003du']


for sn in ['2011fe']:
    # read in data
    data = sncosmo.read_lc(f'data/calibration_sample_lc/{sn}.txt')

    # make model (with dust and redshift)
    # we could probably make the model outside the loop and just update it,
    # but having a clean model brings me piece of mind.
    model = sncosmo.Model(source='salt2', effects=[dust],
                          effect_names=['mw'], effect_frames=['obs'])
    model['mwebv'] = dustmap.ebv(data.meta['ra'], data.meta['dec'])
    model['z'] = data.meta['z']

    # fit data
    print(f'running fitter on {sn}')
    result, fitted_model = sncosmo.mcmc_lc(data, model, ['t0', 'x0', 'x1', 'c'])
    print(f"model parameters: {result.param_names}")
    print(f"best-fit values: {result.parameters}")
    print(f"mean acceptance fraction: {result.mean_acceptance_fraction}")

    # make and save lc
    sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
    plt.savefig(f'figures/lc/{sn}_lc.pdf')
