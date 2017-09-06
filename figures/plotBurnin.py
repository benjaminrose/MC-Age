""" plotBurnin.py
"""
import numpy as np
import matplotlib.pyplot as plt

ID = 'None'
# START = 100
START = 0

# size should be (nwalkers, nsteps, ndim)
samples = np.load('../resources/burnin/SN{}_samples.npy'.format(ID))
# size is (nwalkers, nsteps)
lnprop_resutls = np.load('../resources/burnin/SN{}_lnprob.npy'.format(ID))

print(samples[:,:,0].shape)

f, axarr = plt.subplots(8, sharex=True)
axarr[0].plot(samples[:, START:, 0].T)
axarr[0].set_title('Testing for Burn-in value')
axarr[0].set_ylabel(r'$log_{Z_{sol}}$')
axarr[1].plot(samples[:, START:, 1].T)
axarr[1].set_ylabel(r'$\tau_{dust}$')
axarr[2].plot(samples[:, START:, 2].T)
axarr[2].set_ylabel(r'$\tau$')
axarr[3].plot(samples[:, START:, 3].T)
axarr[3].set_ylabel(r'$t_{0}$')
axarr[4].plot(samples[:, START:, 4].T)
axarr[4].set_ylabel(r'$t_{i}$')
axarr[5].plot(samples[:, START:, 5].T)
axarr[5].set_ylabel(r'$m_{sf}$')
axarr[6].plot(samples[:, START:, 6].T)
axarr[6].set_ylabel('c')
axarr[7].plot(lnprop_resutls.T[START:])
axarr[7].set_ylabel('ln')
# plt.savefig('2017-07-19-burnin.pdf')
plt.show()
