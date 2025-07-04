#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import sys
import os

f=open(sys.argv[1],"r")
s=os.stat(sys.argv[1])

specs=[]
fftlen=n.fromfile(f,dtype=n.uint64,count=1)[0]
# fft length
print(fftlen)
t0=n.fromfile(f,dtype=n.uint64,count=1)[0]
# time of spectrum (unix seconds)
print(t0)
rate=n.fromfile(f,dtype=n.float64,count=1)[0]
print(rate)

f0=n.fromfile(f,dtype=n.float64,count=1)[0]
print(f0)

n_spec=int(n.floor(s.st_size/(8+8+8+8+4*fftlen)-1))
print(n_spec)



S=n.zeros([fftlen,n_spec],dtype=n.float32)
tv=n.zeros(n_spec)
tv[0]=t0
S[:,0]=n.fft.fftshift(n.fromfile(f,dtype=n.float32,count=fftlen))
for i in range(1,n_spec):
    fftlen0=n.fromfile(f,dtype=n.uint64,count=1)[0]
    tv[i]=n.fromfile(f,dtype=n.uint64,count=1)[0]
    rate0=n.fromfile(f,dtype=n.float64,count=1)[0]
    f00=n.fromfile(f,dtype=n.float64,count=1)[0]
    S[:,i]=n.fft.fftshift(n.fromfile(f,dtype=n.float32,count=fftlen))    

max_spec=2000
dt=n.max([1,float(n_spec)/float(max_spec)])
if n_spec < max_spec:
    idx=n.arange(0,n_spec)
else:
    idx=n.array(n.linspace(0,n_spec-1,max_spec),dtype=n.int)

freqs=n.fft.fftshift(n.fft.fftfreq(int(fftlen),d=1.0/float(rate)))+f0
plt.pcolormesh(tv[idx],freqs/1e6,10.0*n.log10(S[:,idx]))
plt.colorbar()
plt.show()


f.close()
