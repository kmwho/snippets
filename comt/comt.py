from __future__ import division
import numpy as np
import pylab as pl
from matplotlib import animation
import nbody

days   = 24*60*60.0  # Assume a normal day for now

def main():
	sys  = nbody.systemFromFile('bodies.dat')
	S    = nbody.RungeKutta4(sys,2433282.5*days) # initial day
	S.integrateTo(2440587.5*days,0.5*days)
	nbody.writeToFile(sys,'output.dat')

if __name__ == '__main__':
	main()
