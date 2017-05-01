""" redoGupta.py -- A file to calculate ages of the whole galaxy, 
like Gupta, with Gupta's galaxies

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-03-23
Python 3.5
"""
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS

def getPhotometry():
    #Get global photometry of objects observed by Gupta
    data = pd.read_csv('data/Gupta11_table2.txt', delimiter='\t', 
                       skiprows=[0,1,2], skipinitialspace=True, 
                       na_values='...', index_col=False)

    #Make a new data table
    HostData = data[['SN ID', 'IAU']]

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
            print(f'host of SN{data.loc[i, "SN ID"]} not found with a spectrum')
            #remove spec option
            results = SDSS.query_region(Galpos, radius=3*u.arcsec, 
                                    photoobj_fields=['objid', 'ra', 'dec', 
                                                     'PetroRad_r', 'u', 'g', 
                                                     'r', 'i', 'z', 'err_u', 
                                                     'err_g', 'err_r', 'err_i',
                                                     'err_z']
                                     )
            if results is None:
                print(f'host of SN{data.loc[i, "SN ID"]} not found')
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

if __name__ == '__main__':
    getPhotometry()