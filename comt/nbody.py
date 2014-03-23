from __future__ import division
import os
import numpy as np
import pylab as pl
from matplotlib import animation

#G = 6.67384e-11 # m3 kg-1 s-2
G = 6.67384e-20 # km3 kg-1 s-2
#G = 1

options = {
	"plot": True
}

class Body:
	def __init__(self, name="Unknown", m = 1.0, r=[.0,.0,.0], v=[.0,.0,.0], a=[.0,.0,.0]):
		self.name = name
		self.r = np.array(r)
		self.v = np.array(v)
		self.a = np.array(a)
		self.m = float(m)
	def cleara(self):
		self.a = np.array([.0,.0,.0])
	def copy(self):
		return Body(self.name, self.m, self.r, self.v, self.a)

class System:
	def __init__(self, bodies=None):
		self.bodies = bodies or []
		self._calcPairs()
	def add(self, bodies):
		self.bodies.extend(bodies)
		self._calcPairs()
	def pairs(self):
		return self._pairs
	def _calcPairs(self):
		self._pairs = [ (self.bodies[i],self.bodies[j]) for (i,j) in Combinations.pairs(len(self.bodies))]
	def copy(self):
		return System( [b.copy() for b in self.bodies] )


class Combinations:
	@staticmethod
	def pairs(n):
		"""Give all possible pairs of numbers from 1 to n
		   TODO change this to iterator"""
		return [(j,i) for i in range(1,n) for j in range(i) ]

class Integrator(object):
	def __init__(self,sys, t=0.0):
		"""Bind instance to System sys"""
		self.sys = sys
		self.t   = t
	def step(self,dt):
		self.t += dt
	def integrateTo(self,t,dt):
		while self.t < t:
			self.step(dt)
		return self.sys

class Symplectic(Integrator):
	def __init__(self, sys, t=0.0):
		"""Symplectic Integrator, extend from this
		and add your order, c and matrices"""
		super(Symplectic,self).__init__(sys, t)

	def step(self,dt):
		"""Step by dt"""
		bodies = self.sys.bodies
		pairs  = self.sys.pairs()
		for i in range(self.order):         #Number of passes is same as the order
			for b in bodies:                #Step Pos
				b.r += b.v*self.c[i]*dt
			[b.cleara() for b in bodies]
			for (b0,b1) in pairs:           #Get predicted accel
				r  = b0.r - b1.r            #TODO seperate derivative functions from integrator
				a  = G/np.dot(r,r)**(3/2)
				b0.a -= a*r*b1.m
				b1.a += a*r*b0.m
			for b in bodies:                #Step Vel
				b.v += b.a*self.d[i]*dt
		self.t += dt

	def integrateTo(self,t,dt):
		super(Symplectic,self).integrateTo(t,dt)
		return self.sys

class Symplectic1(Symplectic):
	def __init__(self, sys, t=0.0):
		"""1st order Symplectic"""
		super(Symplectic1,self).__init__(sys, t)
		self.order = 1
		self.c    = np.array([ 1.0 ])
		self.d    = np.array([ 1.0 ])

class Symplectic4(Symplectic):
	def __init__(self, sys, t=0.0):
		"""4th order Symplectic"""
		super(Symplectic4, self).__init__(sys, t)
		self.order = 4
		k    = 2-(2**(1/3))
		self.c    = np.array([ 1.0/(2*k), (k-1)/(2*k), (k-1)/(2*k), 1.0/(2*k)])
		self.d    = np.array([ 1.0/k, (k-2)/k, 1.0/k, 0.0])

class RungeKutta(Integrator):
	"""Base class for RungeKutta methods
		TODO add all code excepts Butcher Tableau here"""
	def __init__(self,sys, t=0.0):
		super(RungeKutta,self).__init__(sys, t)

class RungeKutta4(RungeKutta):
	"""Data format Y 2d array
	   Each row represents state vector [r,v] of a single body""" 
	def __init__(self,sys, t=0.0):
		super(RungeKutta4,self).__init__(sys, t)
	def _Ydot(self,Y,m):
		# return array of [rdot, vdot]
		Ydot = np.zeros(np.shape(Y))
		pairs      = Combinations.pairs(np.size(Y,0))
		Ydot[:,:3] = Y[:,3:]            # Derivative of r is v
		for (i,j) in pairs:
			r  = Y[i,:3] - Y[j,:3]
			a  = G/np.dot(r,r)**(3/2)
			Ydot[i,3:] -= a*r*m[j]
			Ydot[j,3:] += a*r*m[i]
		return Ydot
	def _Y(self):
		bodies = self.sys.bodies
		Y      = np.empty( (len(bodies), 6) )
		for i in xrange(len(bodies)):
			Y[i,:3] = bodies[i].r
			Y[i,3:] = bodies[i].v
		return Y
	def _updateSystem(self,Y,Ydot=None):
		bodies = self.sys.bodies
		if Ydot is None:
			Ydot = self._Ydot(Y,[b.m for b in bodies])
		for i in xrange(len(bodies)):
			bodies[i].r = Y[i,:3]
			bodies[i].v = Y[i,3:]
			bodies[i].a = Ydot[i,3:]

	def step(self,dt,Y=None,m=None):
		if Y is None:
			Y = self._Y()
			calledDirectly = True
		else:
			calledDirectly = False
		if m is None:
			m = [b.m for b in self.sys.bodies]
		_Ydot = self._Ydot
		k1    = dt*_Ydot(        Y, m)
		k2    = dt*_Ydot( Y+0.5*k1, m)
		k3    = dt*_Ydot( Y+0.5*k2, m)
		k4    = dt*_Ydot(     Y+k3, m)
		Y[:] += k1/6 + k2/3 + k3/3 + k4/6
		self.t += dt
		if calledDirectly:
			self._updateSystem(Y)
		return Y

	def integrateTo(self,t,dt):
		Y     = self._Y()
		m = [b.m for b in self.sys.bodies]
		while self.t < t:
			Y = self.step(dt,Y,m)
		self._updateSystem(Y)
		return self.sys


def systemFromFile(filename):
	bodies = []
	with open(filename,'r') as f:
		for line in f:
			if len(line) == 0 or line[0] == "#":
				continue
			content = line.split(",")
			name    = content[0].strip()
			m       = float(content[1])
			r       = [float(x) for x in content[2:5]]
			v       = [float(x) for x in content[5:8]]
			bodies.append(Body(name,m,r,v))
	return System(bodies)

def writeToFile(system,filename):
	if not os.path.exists('output'):
		os.makedirs('output')
	with open('output/'+filename,'w') as f:
		for b in system.bodies:
			f.write(b.name)
			f.write(', ' + str(b.m))
			f.write(', ' + str(b.r)[1:-1])
			f.write(', ' + str(b.v)[1:-1])
			f.write('\n')

def showAnimation(integ,n=500,dt=1e-3,showEvery=10,**kargs):
	# integ -> Integrator bound to a system
	S   = integ
	sys = S.sys
	tstart = S.t
	fig = pl.figure()
	if 'axlim' in kargs:
		ax  = pl.axes(xlim=kargs['axlim'][0], ylim=kargs['axlim'][1])
	else:
		ax  = pl.axes()
	pl.grid()
	col = kargs['color'] if 'color' in kargs else np.random.random(len(sys.bodies))
	X = [b.r[0] for b in sys.bodies]
	Y = [b.r[2] for b in sys.bodies]
	scat, = pl.plot(X,Y,'o',c=col)
	def update_plot(i,S,scat):
		#dt = 2e-3
		S.integrateTo(i*showEvery*dt+tstart,dt)
		X = [b.r[0] for b in S.sys.bodies]
		Y = [b.r[2] for b in S.sys.bodies]
		scat.set_xdata(X)
		scat.set_ydata(Y)

	ani = animation.FuncAnimation(fig, update_plot, frames=n, interval = 10, fargs=(S,scat))
	ani.save('output/basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
	pl.show()





def main():
	pass


if __name__ == '__main__':
	main()
