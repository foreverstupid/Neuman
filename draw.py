#!/usr/bin/python

from pylab import *

a = loadtxt('graph.plt')


subplot(211)
plot(a[:,0], a[:,1])
plot(a[:,0], a[:,2])

subplot(212)
plot(a[:,0], (a[:,1]-a[:,2]) / a[:,2])

show()
