from __future__ import division
import os
import numpy as np
import pylab as pl
from matplotlib import animation

#G = 6.67384e-11 # m3 kg-1 s-2
G = 1

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
		return [(j,i) for i in range(1,n) for j in range(i) ]

class Integrator(object):
	def __init__(self,sys, t=0.0):
		"""Bind instance to System sys"""
		self.sys = sys
		self.t   = t
	def step(self,dt):
		self.t += dt
	def IntegrateTo(self,t,dt):
		while self.t < t:
			self.step(dt)
		return self.sys

class Symplectic(Integrator):
	def __init__(self, sys, t=0.0):
		"""Symplectic Integrator, extend from this
		and add your order, c and matrices"""
		super(Symplectic,self).__init__(sys, t)

	def step(self,dt):
		bodies = self.sys.bodies
		pairs  = self.sys.pairs()
		for i in range(self.order):
			for b in bodies:
				b.r += b.v*self.c[i]*dt
			[b.cleara() for b in bodies]
			for (b0,b1) in pairs:
				r  = b0.r - b1.r
				a  = G/np.dot(r,r)**(3/2)
				b0.a -= a*r*b1.m
				b1.a += a*r*b0.m
			for b in bodies:
				b.v += b.a*self.d[i]*dt
		self.t += dt

	def IntegrateTo(self,t,dt):
		super(Symplectic,self).IntegrateTo(t,dt)
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
		super(self.__class__, self).__init__(sys, t)
		self.order = 4
		k    = 2-(2**(1/3))
		self.c    = np.array([ 1.0/(2*k), (k-1)/(2*k), (k-1)/(2*k), 1.0/(2*k)])
		self.d    = np.array([ 1.0/k, (k-2)/k, 1.0/k, 0.0])



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
	with open('output'+filename,'w') as f:
		f.write('# System at time '+ sys.t + '\n')
		for b in system.bodies:
			f.write(b.name)
			f.write(', ' + str(b.m))
			f.write(', ' + str(b.r)[1:-1])
			f.write(', ' + str(b.v)[1:-1])
			f.write('\n')


### Everything below is hardcoded

def plotSystem(sys,**kargs):
	if not options["plot"]:
		return
	X = [b.r[0] for b in sys.bodies]
	Y = [b.r[1] for b in sys.bodies]
	return pl.scatter(X,Y,**kargs)

def showPlot():
	if not options["plot"]:
		return
	pl.grid()
	pl.show()

def showAnimation():
	sys = systemFromFile('bodies.dat')
	S   = Symplectic4(sys)
	fig = pl.figure()
	pl.grid()
	ax  = pl.axes(xlim=(-1.5, 1.5), ylim=(-0.5, 0.5))
	col = np.random.random(3)
	X = [b.r[0] for b in sys.bodies]
	Y = [b.r[1] for b in sys.bodies]
	scat, = pl.plot(X,Y,'o',c=col)
	def update_plot(i,S,scat):
		dt = 2e-3
		S.IntegrateTo(i*10*dt,dt)
		X = [b.r[0] for b in S.sys.bodies]
		Y = [b.r[1] for b in S.sys.bodies]
		scat.set_xdata(X)
		scat.set_ydata(Y)

	ani = animation.FuncAnimation(fig, update_plot, frames=200, interval = 10, fargs=(S,scat))
	ani.save('basic_animation.mp4', fps=12, extra_args=['-vcodec', 'libx264'])
	pl.show()


def main():
	sys = systemFromFile('bodies.dat')
	S   = Symplectic4(sys)
	plotSystem(sys,s=80,marker='o')
	for i in xrange(1000):
		S.IntegrateTo(i,1e-1)
		plotSystem(sys,s=80,marker='+')
	showPlot()
	writeToFile(sys,'output/output.dat')


if __name__ == '__main__':
	showAnimation()
