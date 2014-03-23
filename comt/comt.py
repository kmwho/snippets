from __future__ import division
import numpy as np
from numpy.linalg import norm
import pylab as pl
from matplotlib import animation
import nbody, plot3D

days   = 24*60*60.0  # Assume a normal day for now

class RK4adaptive(nbody.RungeKutta4): # Currently only checking for the probe
	def integrateTo(self,t,dt):
		Y     = self._Y()
		m = [b.m for b in self.sys.bodies]
		Jupiter = self.sys.bodies[1]
		probe   = self.sys.bodies[2]   # The Probe
		self._updateSystem(Y,self._Ydot(Y,m))
		while self.t < t:
			delt = dt
			#while True:
				#vn     = probe.a*delt          # Simply checking if the velocity is deviating a lot
				#check  = norm(vn)/norm(probe.v)# It is much better to use jerk, snap etc., but that takes some work
				#if check < 0.01:
				#	break
				#else:
				#	print check
				#	delt /= 2
			if norm(Jupiter.r - probe.r) < 1e6:
				delt /= 10
				print "in Jupiter SOI"
			if np.abs(delt-dt) > 0.1*dt:
				print dt/delt
			Y = self.step(delt,Y,m)
			self._updateSystem(Y)
		return self.sys

def main():
	sys  = nbody.systemFromFile('bodies.dat')
	S    = nbody.RungeKutta4(sys,2433282.5*days) # initial day
	S.integrateTo(2440587.5*days,0.5*days)
	#nbody.showAnimation(S,7305//5,1*days,5,axlim=[(-1e9,1e9),(-1e9, 1e9)]) #Simulation is for 7305 days
	nbody.writeToFile(sys,'output.dat')

def main2():
	sys  = nbody.systemFromFile('bodies.dat')
	S    = RK4adaptive(sys,2433282.5*days) # initial day
	#S.integrateTo(2440587.5*days,0.5*days)
	#nbody.showAnimation(S,7305//10,2*days,5,axlim=[(-1e9,1e9),(-1e9, 1e9)]) #Simulation is for 7305 days
	plot3D.showAnimation(S,7305//10,2*days,5,axlim=[(-1e9,1e9),(-1e9, 1e9),(-1e9, 1e9)])
	nbody.writeToFile(sys,'output.dat')


if __name__ == '__main__':
	main2()
