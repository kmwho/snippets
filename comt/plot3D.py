import pylab as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

def showAnimation(integ,n=500,dt=1e-3,showEvery=10,**kargs):
	# integ -> Integrator bound to a system
	S   = integ
	sys = S.sys
	tstart = S.t
	fig = pl.figure(1)
	fig.clf()
	ax = fig.add_subplot(111, projection='3d')
	if 'axlim' in kargs:
		ax = fig.add_subplot(111, projection='3d',xlim=kargs['axlim'][0], ylim=kargs['axlim'][1], zlim=kargs['axlim'][1])
	else:
		ax = fig.add_subplot(111, projection='3d')
	#ax.grid(True,'both')
	col = kargs['color'] if 'color' in kargs else np.random.random(len(sys.bodies))
	X = [b.r[0] for b in sys.bodies]
	Y = [b.r[1] for b in sys.bodies]
	Z = [b.r[2] for b in sys.bodies]
	scat, = pl.plot(X,Y,Z,'o',c=col)
	def update_plot(i,S,scat):
		#dt = 2e-3
		S.integrateTo(i*showEvery*dt+tstart,dt)
		X = [b.r[0] for b in S.sys.bodies]
		Y = [b.r[1] for b in S.sys.bodies]
		Z = [b.r[2] for b in S.sys.bodies]
		scat.set_xdata(X)
		scat.set_ydata(Y)
		scat.set_3d_properties(Z)

	ani = animation.FuncAnimation(fig, update_plot, frames=n, interval = 10, fargs=(S,scat))
	ani.save('output/basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
	pl.show()

