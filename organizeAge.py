""" organizeAge.py -- this grabs MCMC results and aggregates it into one tsv

Benjamin Rose
brose3@nd.edu/benjamin.rose@me.com
Notre Dame
Python 3.5
2017-05-23
"""

import glob
import re

import numpy as np
from astropy.table import Table
import pandas as pd

#get files you want to make tsv from
fileNames = glob.glob('resources/global/SN*.tsv')

#make collective tsv holding DataFrame
# tableData = pd.DataFrame(columns=['SN', 'age', 'age-uncert', 'age+uncert'])
tableData = Table(names=('SN', 'age', 'age-uncert', 'age+uncert'), dtype=('i', 'f', 'f', 'f'))

#set up Regular Expression (for getting SN number) once
p = re.compile(r'\d+')

#for each file
i = 1
for f in fileNames:
    print('Working on {}'.format(i))
    #import mcmc data
    data = Table.read(f, format='ascii.commented_header', delimiter='\t')
    data = data.to_pandas()

    #calculate age percentiles
    age_precentiels = list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), 
                               np.nanpercentile(data['age'], [16, 50, 84]).reshape(1,3)))[0]

    #get SN number
    number = p.findall(f)     #should return a list of 1 string

    #save to holding DataFrame
    #convert list of a string to list of an int
    hold = [int(number[0])]
    hold.extend(age_precentiels)

    tableData.add_row(hold)

    i += 1
    #debugging breakout
    # if i > 5:
    #     break

#save DataFrame to a tsv
print('\n')
tableData['age'].unit = 'Gyr'
tableData.sort('SN')
print(tableData)

tableData.write('del.tsv', format='ascii.tab', overwrite=True)