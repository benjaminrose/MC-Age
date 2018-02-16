""" redoGupta.py -- A file to calculate ages of the whole galaxy, 
like Gupta, with Gupta's galaxies

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-03-23
Python 3.5
"""
import logging

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

import calculateAge

module_logger = logging.getLogger("fsps-age.redoGupta")

def getPhotometry(dataset='Gupta'):
    """use astroquery to get host galaxy global photometry 

    Parameters
    ----------
    dataset: str
        Do you want to get the data for calibration (`"Gupta"`) or to compare
        local to global values (`"Campbell"`).
    """
    #this gives annoying warnings, so I moved it here
    from astroquery.sdss import SDSS

    if dataset == "Gupta":
        #Get objects analyzed by Gupta
        data = pd.read_csv('data/Gupta11_table2.txt', delimiter='\t', 
                           skiprows=[0,1,2], skipinitialspace=True, 
                           na_values='...', index_col=False)
    elif dataset == "Campbell":
        #Get Campbell/Holtzmann objects & set up to match Gupta file
        data = pd.read_csv('data/CampbellHoltzman.tsv', delimiter='\t')
        data = data.rename(columns={'SNID': 'SN ID'})
        data['RA'] = np.nan
        data['Dec'] = np.nan

        # Get data set with RA & Dec from Campbell's data set
        hdu = fits.open('data/SDSS_Photometric_SNe_Ia.fits')
        campbell = hdu[1].data

        # Get RA and Dec from other file
        for i in range(len(data)):
            # get index of campbell for next SN in data
            index = np.where(campbell['CID'] == data.loc[i, 'SN ID'])[0]

            #store RA & dec
            data.loc[i, 'RA'] = campbell['GAL_RA'][index]
            data.loc[i, 'Dec'] = campbell['GAL_DEC'][index]

    else:
        #todo(convert to error, log, or warning)
        print("can only look up 'Gupta' or 'Campbell' data sets")
        from sys import exit; exit()

    #Make a new data table
    # HostData = data[['SN ID', 'IAU']]

    for i in range(len(data)):
        Galpos = SkyCoord(ra=data.loc[i, 'RA']*u.deg, 
                         dec=data.loc[i, 'Dec']*u.deg)
        # Galpos = SkyCoord(ra=data.loc[0, 'RA']*u.deg, 
                         # dec=data.loc[0, 'Dec']*u.deg)
        print(data.loc[i, 'SN ID'])
        radius=5*u.arcsec
        results = SDSS.query_region(Galpos, spectro=True, radius=radius, 
                                    specobj_fields=['z'], 
                                    photoobj_fields=['objid', 'ra', 'dec', 
                                                     'PetroRad_r', 'u', 'g', 
                                                     'r', 'i', 'z', 'err_u', 
                                                     'err_g', 'err_r', 'err_i',
                                                     'err_z']
                                     )
        if results is None:
            print('host of SN{} not found with a spectrum'.format(data.loc[i, "SN ID"]))
            #remove spec option
            results = SDSS.query_region(Galpos, radius=3*u.arcsec, 
                                    photoobj_fields=['objid', 'ra', 'dec', 
                                                     'PetroRad_r', 'u', 'g', 
                                                     'r', 'i', 'z', 'err_u', 
                                                     'err_g', 'err_r', 'err_i',
                                                     'err_z']
                                     )
            if results is None:
                print('host of SN{} not found'.format(data.loc[i, "SN ID"]))
                continue
        # import pdb; pdb.set_trace()
        if len(results) == 1:
            # check redshift
            # just save the one result to HostData
            #print the whole Table, no matter how many characters wide
            results.pprint(max_width=-1)
        if len(results) > 1:
            firstID = results[0]['objid']
            isSame = True
            for j in results:
                if firstID != j['objid']:
                    isSame = False
            if isSame:
                print(results[0])#.table.pprint(max_width=-1)
            else:
                print('length:', len(results))
                results.pprint(max_width=-1)


#Run `calculateAge.py` running 15 objects (~5 days on CRC) per job array 
#Gupta has ~210 objects so a complete job array needs 14 objects

# todo this might need a new name.
def redoGupta(jobID, lenJobs=50, debug=False, dataset='circle'):
    """This runs through each broken-up tsv of the global photometry and calling
    calculateAge to save global age calculations.
    
    Parameters
    ----------
    jobID : int
        This is the ID of the global photometry tsv to be used. Should come from
        the crc job-array ID. CRC says it NEEDS to be 1 indexed.

    lenJobs : int
        The total number of jobs that will be run when "embarrassingly"
        parallelize through the input file. This lets each `jobID` know what
        line of data it should analyze. This is not the max value of jobID, 
        but the total number of jobID's, like `len()`.

    debug : bool
        Flag to have MCMC run incredibly short and in no way accurately. Also 
        does not save resulting chain. Should take around ~12 mins to get a 
        value of one SN.

    dataset : str, ['gupta', 'messier', 'circle', 'campbell']
        flag to switch between two "calibration" methods: redoing Gupta's
        analysis (with `'gupta'`), looking at about ~10 local galaxies and making sure spirals are younger then ellipticals (with `'messier'`), or testing to recover "known" star formation histories (with `'circle'). This defaults to `'circle'`. More details on circle photometry can be found at . . .
    """
    logger = logging.getLogger("fsps-age.redoGupta.redoGupta")

    # Import data file
    if dataset in ['gupta', 'messier', 'circle', 'campbell', 'campbellG']:
        if dataset == 'gupta':
            logger.info('importing GlobalPhotometry-Gupta.tsv')
            data = pd.read_csv('data/GlobalPhotometry-Gupta.tsv',
                                delimiter='\t', skiprows=[0,1,2,4],
                                skipinitialspace=True, na_values='...',
                                index_col=False)
        elif dataset == "messier":
            # these galaxies are more "known' then gupta, at least by type
            logger.info('importing galaxyPhotometry.tsv')
            data = pd.read_csv('data/galaxyPhotometry.tsv', delimiter='\t')
        elif dataset == 'campbell':
            logger.info('importing CampbellHoltzman.tsv')
            data = pd.read_csv('data/CampbellHoltzman.tsv', delimiter='\t')
        elif dataset == 'campbellG':
            logger.info('importing CampbellHoltzman_global.tsv')
            data = pd.read_csv('data/CampbellHoltzman_global.tsv',
                               delimiter='\t')
        else:
            # default to circle test
            logger.info('importing circlePhotometry.tsv')
            data = pd.read_csv('data/circlePhotometry.tsv', delimiter='\t')
    else:
        raise ValueError("Invalid dataset argument")
    
    #cut down dataset
    Ndata = len(data)
    stepSize = int(np.ceil(Ndata/lenJobs))
    startIndex = (int(jobID)-1)*stepSize
    endIndex = startIndex + stepSize

    #check if start value is outside data size
    if startIndex >= Ndata:
        logger.info('This jobID is not needed. Start value larger than data' +
                    ' size')
        return  

    #use `.loc` or else zips won't work
    if endIndex >= Ndata:
        #using `[x:x]` returns a DataFrame, Using `[x]` returns a Series
        data = data.loc[startIndex:Ndata-1]
    else:
        #`loc` is inclusive, WHY!!
        data = data.loc[startIndex:endIndex-1]

    # Iterate over each SN in data file
    #iterate over a zip of SNID, photometry (zipped together), photometry 
    #uncertainty (zipped together), and redshift
    #note http://stackoverflow.com/questions/10729210/iterating-row-by-row-through-a-pandas-dataframe 
    #and http://stackoverflow.com/questions/7837722/what-is-the-most-efficient-way-to-loop-through-dataframes-with-pandas/34311080#34311080
    for sn, photometry, uncertainty, redshift in zip(
                data['SNID'], 
                zip(data['u'], data['g'], data['r'], data['i'], data['z']),
                zip(data['err_u'], data['err_g'], data['err_r'], data['err_i'],
                    data['err_z']),
                data['redshift']
            ):

        # calculate age
        #does this save the data tagged as global?
        logger.info('getting age for SN' + str(sn))
        age = calculateAge.calculateAge(redshift, photometry, uncertainty, 
                                        SNID=sn, debug=debug, dataset=dataset)
        logger.info("Global age for SN{}: ".format(sn) + str(age))
        print("Global age for SN{}: ".format(sn), age)
        

if __name__ == '__main__':
    getPhotometry('Campbell')
    # redoGupta(1)