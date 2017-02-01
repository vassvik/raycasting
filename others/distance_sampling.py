import numpy
import matplotlib.pyplot as plt 



def getBilinear(V00, V10, V01, V11, Nx, Ny, f):
	mat = numpy.zeros((Ny, Nx))

	for j in range(Ny):
		for i in range(Nx):
			x = f(i/(Nx-1.0))
			y = f(j/(Ny-1.0))
			mat[j][i] = V00*(1.0-x)*(1.0-y) + V10*x*(1.0-y) + V01*(1.0-x)*y + V11*x*y
	return mat

N = 128
cmap = 'nipy_spectral'
cmap = 'Paired'



funs = [
	lambda t: numpy.floor(1.9999*t), 
	lambda t: t,
	lambda t: t**2.0/(t**2.0 + (1-t)**2.0),
	lambda t: t**3.0/(t**3.0 + (1-t)**3.0),
	lambda t: t**6.0/(t**6.0 + (1-t)**6.0),
	lambda t: -2*t**3 + 3*t**2,
	lambda t: 10*t**3 - 15*t**4 + 6*t**5,
	lambda t: 0.5 - 0.5*numpy.cos(numpy.pi*t)
]
strings = [
	r"$H(x-1/2)$",
	r"$t$",
	r"$t^2/(t^2 + (1-t)^2)$",
	r"$t^3/(t^3 + (1-t)^3)$",
	r"$t^6/(t^6 + (1-t)^6)$",
	r"$3t^2 - 2 t^3$",
	r"$10 t^3 - 15 t^4 + 6 t^5$",
	r"$(1- \cos(\pi t))/2$"

]


size = (21,12)
"""
Nx = 10
Ny = 10
biggrid = numpy.random.randint(0, 2, size=(Nx, Ny))
for i in range(Nx):
	for j in range(Ny):
		if i == 0 or j == 0 or i == Nx-1 or j == Ny-1:
			biggrid[j][i] = 0

plt.figure(figsize=size)
for k in range(len(funs)):
	plt.subplot(2,4,k+1)
	grid = numpy.zeros((N*(Nx-1), N*(Ny-1)))
	for i in range(Nx-1):
		for j in range(Ny-1):
			grid[j*N:(j+1)*N,(i+0)*N:(i+1)*N] = getBilinear(biggrid[j,i], biggrid[j,i+1], biggrid[j+1,i], biggrid[j+1,i+1], N, N, funs[k])
	plt.imshow(grid, interpolation="None")
	plt.title(strings[k])
	plt.colorbar()
plt.savefig("v.png")
"""

a = 1.0

plt.figure(figsize=size)
x = [1.0*i/999 for i in range(1000)]
t = numpy.power(x,a)
for k in range(len(funs)):
	print(k)
	y = funs[k](t)
	plt.plot(x, y, label=strings[k])
plt.legend(loc='best')
plt.savefig("f.png")

plt.figure(figsize=size)
for k in range(len(funs)):
	y = funs[k](t)
	dy = numpy.gradient(y)
	dy = dy/max(dy)
	plt.plot(x, dy, label=strings[k])
plt.legend(loc='best')
plt.savefig("dx.png")

plt.figure(figsize=size)
for k in range(len(funs)):
	y = funs[k](t)
	if strings[k] == r"$t$":
		dy2 = 0*numpy.gradient(numpy.gradient(y))
	else:
		dy2 =  numpy.gradient(numpy.gradient(y))
		dy2 = dy2 / max(dy2)
	plt.plot(x, dy2, label=strings[k])
plt.legend(loc='best')
plt.savefig("dx2.png")

plt.show(block=True)