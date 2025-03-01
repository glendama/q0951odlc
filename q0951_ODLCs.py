#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

t, A, Aerr, B, Berr = np.loadtxt('q0951LT_USNO_DES_PS_2008_2023.dat', unpack=True)
N = len(t)
print (N)

dt1 = 16.0	# time delay for approach 1 
dt2 = 13.3	# time delay for approach 2 
delta = 10	# time window for selection/rejection 
orderfit = 1	# degree of the fitting polynomial 

##############
# Approach 1 #
##############
 
# INTERPOLATION OF THE SHIFTED LIGHT CURVE OF A  
Aint1 = np.interp(t, t+dt1, A)
Aint1_err = np.interp(t, t+dt1, Aerr)
BAint1 = B - Aint1
BAint1_err = np.sqrt(Berr**2 + Aint1_err**2)

# DATA SELECTION & ODLC
ind = np.zeros(len(t))
for i in range(N-1):
    for j in range(N):
        if (t[i]+dt1-t[j]) < 0 and (t[i+1]+dt1-t[j]) > 0 and ((t[i]+dt1-t[j]) > -delta or (t[i+1]+dt1-t[j]) < delta):
            ind[j] = 1
tBA1 = np.compress(ind, t)            
BA1 = np.compress(ind, BAint1)
BA1err = np.compress(ind, BAint1_err)

BA1ave = np.average(BA1)
np.savetxt('q0951_ODLC1.dat', np.c_[tBA1, BA1-BA1ave, BA1err], fmt='%.3f') 

# FITTING A POLYNOMIAL MODEL TO THE ODLC
z = np.polyfit(tBA1, BA1-BA1ave, orderfit)
p = np.poly1d(z)
BA1fit = p(tBA1)

##############
# Approach 2 #
##############
 
# INTERPOLATION OF THE SHIFTED LIGHT CURVE OF A  
Aint2 = np.interp(t, t+dt2, A)
Aint2_err = np.interp(t, t+dt2, Aerr)
BAint2 = B - Aint2
BAint2_err = np.sqrt(Berr**2 + Aint2_err**2)

# DATA SELECTION & ODLC
ind = np.zeros(len(t))
for i in range(N-1):
    for j in range(N):
        if (t[i]+dt2-t[j]) < 0 and (t[i+1]+dt2-t[j]) > 0 and ((t[i]+dt2-t[j]) > -delta or (t[i+1]+dt2-t[j]) < delta):
            ind[j] = 1
tBA2 = np.compress(ind, t)            
BA2 = np.compress(ind, BAint2)
BA2err = np.compress(ind, BAint2_err)

BA2ave = np.average(BA2)
np.savetxt('q0951_ODLC2.dat', np.c_[tBA2, BA2-BA2ave, BA2err], fmt='%.3f') 

# FITTING A POLYNOMIAL MODEL TO THE ODLC
z = np.polyfit(tBA2, BA2-BA2ave, orderfit)
p = np.poly1d(z)
BA2fit = p(tBA2)

tmin = 54466; tmax = 60309; t0 = 54466; year0 = 2008; Nyears = 16
ymin, ymax = 17.90, 17.12
years = year0 + np.arange(Nyears)
years_firstday = t0 + 365.25*np.arange(Nyears)
ms = 5
mec = 'k'; mew = 1

# Figure 1: LIGHT CURVES 
plt.errorbar(t, A, yerr=Aerr, fmt='o', color='r', mew=mew, mec=mec, ms=ms, label=r'q0951A')
plt.errorbar(t, B-1.1, yerr=Berr, fmt='s', color='b', mew=mew, mec=mec, ms=ms, label=r'q0951B-1.1')
plt.legend()
plt.xlabel(r'MJD (days)', size='xx-large')
plt.ylabel(r'$r$ (mag)', size='xx-large')
plt.xlim(tmin, tmax)
plt.ylim(ymin, ymax)
ax2 = plt.twiny()
ax2.set_xlim(tmin, tmax)
ax2.set_xticks(years_firstday)
ax2.set_xticklabels((years), ha='left',fontsize=7)
ax2.set_xlabel(r'year', fontsize='xx-large')
plt.savefig('q0951_LCs.png', dpi=300)
plt.show()

# Figure 2: ODLCs & linear fits
plt.errorbar(tBA1, BA1-BA1ave, yerr=BA1err, fmt='o', color='g', ms = 7, mec=mec, mew=mew, label='ODLC1')
plt.errorbar(tBA2, BA2-BA2ave, yerr=BA2err, fmt='s', color='r', ms = 5, mec=mec, mew=mew, label='ODLC2')
plt.legend()
plt.plot(tBA1, BA1fit, 'C2--')
plt.plot(tBA2, BA2fit, 'C3--')
plt.xlabel(r'MJD (days)', size='xx-large')
plt.ylabel(r'ODLC (mag)', size='xx-large')
plt.xlim(tmin, tmax)
ymin, ymax = 0.25, -0.25
plt.ylim(ymin, ymax)
ax2 = plt.twiny()
ax2.set_xlim(tmin, tmax)
ax2.set_xticks(years_firstday)
ax2.set_xticklabels((years), ha='left',fontsize=7)
ax2.set_xlabel(r'year', fontsize='xx-large')
plt.savefig('q0951_ODLCs.png', dpi=300)
plt.show()
