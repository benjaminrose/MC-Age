# Usage

## burnin

## run

currently only runs global tests (on SN hosts analsyied by Gupta 2011 or some Messier galaxies). Also runs a circle test.

### Circle

The goal of this is to test the MCMC run. We can put in some star formation parameters into FSPS and get out some magnitudes, can `fsps-age` properly use MCMC to retrieve these star formation parameters?

Using the stellar popultion [used in the analsysis](https://github.com/benjaminrose/SNIa-Local-Environments/blob/3ee9b1f508c2583bcd8a314827ab7c0497526f93/calculateAge.py#L293), I was able to get a set of SED's. The resulting "data" can be found as [`circlePhotometry.tsv`](https://github.com/benjaminrose/SNIa-Local-Environments/blob/master/data/circlePhotometry.tsv). The details of its creation are below.

We used the operation in [the code](https://github.com/benjaminrose/SNIa-Local-Environments/blob/3ee9b1f508c2583bcd8a314827ab7c0497526f93/calculateAge.py#L84) of `dust1 = 2.0*dust2`.

We choose our redshift to be 0.05 resulting in the age of emission being 13.18 Gyr using the [standard cosmology](https://github.com/benjaminrose/SNIa-Local-Environments/blob/3ee9b1f508c2583bcd8a314827ab7c0497526f93/calculateAge.py#L43).

ID |  comment         | log(z_sol) | tau_dust | tau | t_start | t_transition | m_sf
---|------------------|------------|----------|-----|---------|--------------|-------
1  | old w/ metals    |    -0.5    |   0.1    | 0.5 |   1.5   |      9.0     | -1.0
2  | old & young      |    -0.5    |   0.1    | 0.5 |   1.5   |      9.0     | 15.0
3  | young & younger  |    -0.5    |   0.1    | 7.0 |   3.0   |      10      | 15.0
4  | young            |    -0.5    |   0.1    | 7.0 |   3.0   |      13      | 0
5  | old & metal poor |    -1.5    |   0.1    | 0.5 |   1.5   |      9.0     | -1.0
6  | young & dusty    |    -0.5    |   0.8    | 7.0 |   3.0   |      10      | 15.0

ID |  comment         | u - g  |   u   |   g   |   r   |   i   |  z
---|------------------|--------|-------|-------|-------|-------|-------
1  | old w/ metals    |  1.60  | 45.36 | 43.76 | 42.99 | 42.67 | 42.39
2  | old & young      |  1.57  | 45.31 | 43.74 | 42.98 | 42.66 | 42.39
3  | young & younger  |  0.72  | 41.15 | 40.43 | 40.40 | 40.19 | 40.21
4  | young            |  0.91  | 42.65 | 41.74 | 41.49 | 41.26 | 41.16
5  | old & metal poor |  1.40  | 44.69 | 43.29 | 42.70 | 42.45 | 42.29
6  | young & dusty    |  1.08  | 42.66 | 41.58 | 41.25 | 41.01 | 40.86

Adding a c = -25 we get the SED's in [`circlePhotometry.tsv`](https://github.com/benjaminrose/SNIa-Local-Environments/blob/master/data/circlePhotometry.tsv).