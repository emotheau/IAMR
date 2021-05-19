import numpy as np
import csv
import matplotlib.pyplot as plt
import os

basedir = '/project/projectdirs/dasrepo/jpathak/IAMR-HIT/Exec/HIT/ANALYSIS/'
#plt_00475-slice-spectra.csv
filenames = ['plt_00000-slice-spectra.csv', 'plt_00393-slice-spectra.csv', 'plt_00475-slice-spectra.csv']
ii = 0
for filename in filenames:

    plt.figure()
    with open(os.path.join(basedir, filename)) as csvfile:
        specreader = csv.reader(csvfile)
        for row in specreader:
            print(row)
    
    spectrum1 = np.genfromtxt(os.path.join(basedir, filename), delimiter=',')
    ii+=1
    plt.loglog(spectrum1[1:,0], spectrum1[1:, 1], basex=2, basey=2)
    plt.xlabel("K")
    plt.ylabel("p(K)")
    plt.savefig('spectrum_128_' + str(ii))

