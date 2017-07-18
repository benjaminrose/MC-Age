""" plotBurnin.py
"""
import numpy as np
import matplotlib.pyplot as plt

# size should be (nwalkers, nsteps, ndim)
samples = np.load('../resources/burnin/samples.npy')
# size == (nwalkers, nsteps)
# lnprop_resutls = np.load('../resources/burnin/lnprob.npy')


f, axarr = plt.subplots(8, sharex=True)
axarr[0].plot(samples[:,:,0].T)
axarr[0].set_title('Testing for Burn-in value')
axarr[0].set_ylabel('logzsol')
axarr[1].plot(samples[:,:,1].T)
axarr[1].set_ylabel('dust')
axarr[2].plot(samples[:,:,2].T)
axarr[2].set_ylabel('tau')
axarr[3].plot(samples[:,:,3].T)
axarr[3].set_ylabel('tStart')
axarr[4].plot(samples[:,:,4].T)
axarr[4].set_ylabel('sfTrans')
axarr[5].plot(samples[:,:,5].T)
axarr[5].set_ylabel('sfSlope')
axarr[6].plot(samples[:,:,6].T)
axarr[6].set_ylabel('c')
# axarr[7].plot(lnprop_resutls.T)
# axarr[7].set_ylabel('ln')
# plt.savefig('figures/burnin.pdf')
plt.show()
##################