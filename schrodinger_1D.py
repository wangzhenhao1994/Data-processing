#!/usr/bin/env python
import matplotlib.pyplot as plt

from classes.class_secant import Secant_Method
from classes.class_numerov import Numerov_Method
from classes.class_simpson import Simpson_Method


'''
Potential function of the quantum system.
---------------------------------------------------------
'''

from numpy import cosh
def v(x):
	'This is the function where analytical form is given'
	h_bar = 1.0
	m = 1.0
	alpha = 1.0
	lamb = 4.0

	frst = (h_bar/(2*m))*alpha*alpha * lamb*(lamb - 1)
	scnd = 0.5 - 1/(cosh(alpha*x)*cosh(alpha*x))

	return frst * scnd
'''
def v(x):
	'Square Well Potential'
	if abs(x) > 2:
		return 0
	else:
		return x*x - 4
'''

'''
Help function to discretization of the analytical form
---------------------------------------------------------
'''

def anal_disc():
	'Discretization of the analytical form will be done'
	initial = -10
	h = 0.01

	x = []
	y = []

	while initial <= 10:
		x.append(initial)
		y.append(v(initial))
		initial += h

	return x, y

'''
Function used to normalized the wavefunction.
---------------------------------------------------------
'''

def normalized(u):
	'Normalization to the wavefunction'

	result = integration.simpson(u)

	for i in range(len(u)):
		u[i] = u[i] / result

	return u

'''
Function used to make list of the eigenenery to plot.
---------------------------------------------------------
'''

def eigen(x, l):
	'Making the list of an eigenvalue for plotting'
	y = []
	for i in range(len(x)):
		y.append(l)

	return y

'''
Function used to add eigenfunction with eigenenergy.
---------------------------------------------------------
'''

def mixfunction(u, el):
	'For nice visualization adding eigenfunction with eigenvalue'
	y = []
	for i in range(len(u)):
		y.append(u[i]+el[i])
	return y


'''
Stepping into the main program.
---------------------------------------------------------
'''

'Plotting potential'
x, v = anal_disc()
plt.plot(x,v, 'k')

'Calling Numerov Method for the potential v(x)'
wavefunction = Numerov_Method(x, v)

'Calling Secant method for the function $numerov$'
find = Secant_Method(wavefunction.numerov)

'Calling Simpson Method to integrate'
integration = Simpson_Method()

'Initial guess for eigen value'
l = -10

title = 'Eigenvalue(s) = '

#plt.figure(1)

while l < max(v):
	'Searching for the eigenvalue.'
	l = find.secant(l)

	title = title + str("%.2f" % round(l,2)) + ', '

	'Calculating the eigenfunction.'
	x, u = wavefunction.numerov(l)

	'Normalizing the eigenfunction'
	u = normalized(u)

	'list to plot eigenvalue'
	el = eigen(x, l)

	'Adding eigenvalue and normalized wavefunction'
	n_u = mixfunction(u, el)

	plt.plot(x, n_u, 'g')
	plt.plot(x, el, 'r')

	'for loop'
	l += 0.3

'Plotting graph'
plt.figure(1)
plt.xlabel('Position x')
plt.ylabel('Magnitude')
plt.title(title)
plt.legend(['$V(x)$', '$\phi(x)$', '$E_n$'], loc=4)
plt.axis([-5,5,min(v)-0.5,max(v)+0.5])

plt.show()