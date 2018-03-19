#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 11:24:54 2018

@author: chenhao
"""
import matplotlib.pyplot as plt
import numpy as np

sarprsity = np.array([	24.13,	23.54,	24.04, 26.1996, 34.65])
lrl2 =	np.array([0.9913,	1,		0.9993,	1, 		0.9927])
lrl3 =	np.array([0.9903,	0.9933,	0.9893, 0.9950, 0.9827])
ZFl2 =	np.array([0.964,	0.9393,	0.9747,	0.9553, 0.3147])
ZFl3 =	np.array([0.8287,	0.8283,	0.8543, 0.8437, 0.4857])


plt.close('all')
 
plt.figure()
plt.plot(sarprsity)
plt.plot(sarprsity, 'x')
plt.minorticks_on()
plt.grid(b = True, which='minor',linestyle='--')
plt.title('sparsity level for different dictionaries')
plt.ylabel('Sparsity level')
plt.xticks(range(4), ('FDDLO', 'beta=0.1', 'beta=1', 'beta=5', 'beta=10'))

plt.figure()
plt.plot(lrl2)
plt.plot(lrl3)
plt.plot(ZFl2)
plt.plot(ZFl3)
plt.plot(lrl2, 'x')
plt.plot(lrl3, 'x')
plt.plot(ZFl2, 'x')
plt.plot(ZFl3, 'x')
plt.minorticks_on()
plt.grid(b = True, which='both',linestyle='--')
plt.title('classification accuracy for different dictionaries')
plt.legend(['logistic regression for L =2','logistic regression for L =3', 'Zero forcing for L =2', 'Zero forcing for L =3'])
plt.ylabel('accuracy')
plt.xticks(range(5), ('FDDLO', 'beta=0.1', 'beta=1', 'beta=5','beta=10'))

plt.show()