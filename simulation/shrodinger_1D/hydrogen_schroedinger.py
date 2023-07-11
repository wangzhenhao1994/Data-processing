#!/usr/bin/env python
import matplotlib.pyplot as plt

from classes.class_secant import Secant_Method
from classes.class_numerov import Numerov_Method
from classes.class_simpson import Simpson_Method

'''
Please provide the input value for solving SE of Hydrogen
Atom
=========================================================
eigen_min is the lowest probable eigen value
potential_initial_position is the lowest position of potential
potential_final_position is the highest position of potential
number_of_eigen is the number of eigenfunction
=========================================================
'''

angular_l = int(input("Please provide the total angular momentum: "))

number_of_eigen = 5
eigen_min = -10
potential_initial_position = 0.01
potential_final_position = 20.0

'''
Potential function of the quantum system.
---------------------------------------------------------
'''


def v(x):
    'potential for hydrogen atom'

    return -1.0/x + angular_l * (angular_l + 1.0)/(2*x*x)


'''
Help function to discretization of the analytical form
---------------------------------------------------------
'''


def numeric_potential():
    'Discretization of the analytical form will be done'
    initial = potential_initial_position
    h = 0.01

    x = []
    y = []

    while initial <= potential_final_position:
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
Function used to calculate the probability
--------------------------------------------------------
'''


def probability(u):
    y = []
    for i in range(len(u)):
        y.append(u[i] * u[i])

    result = integration.simpson(y)

    for i in range(len(u)):
        y[i] = y[i]/result

    return y


'''
Function used to make list of the eigenenery to plot.
---------------------------------------------------------
'''


def eigen(x, eigen_value):
    'Making the list of an eigenvalue for plotting'
    y = []
    for i in range(len(x)):
        y.append(eigen_value)

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
x, v_x = numeric_potential()
plt.plot(x, v_x, 'k')


'Calling Numerov Method for the potential v(x)'
wavefunction = Numerov_Method(x, v_x)

'Calling Secant method for the function $numerov$'
find = Secant_Method(wavefunction.numerov)

'Calling Simpson Method to integrate'
integration = Simpson_Method()

'Figure title'
title = 'Eigenvalue(s) = '

y_min = 0
y_max = 0
count = 0
while (count < number_of_eigen and eigen_min < 0):
    'Searching for the eigenvalue.'
    eigen_value = find.secant(eigen_min)

    if(eigen_value < 0):
        title = title + str("%.3f" % round(eigen_value, 3)) + ', '

        'Calculating the eigenfunction.'
        x, u = wavefunction.numerov(eigen_value)

        'Normalizing the eigenfunction'
        u = normalized(u)

        'Probability of eigenfunction'
        p = probability(u)

        'list to plot eigenvalue'
        el = eigen(x, eigen_value)

        'Adding eigenvalue and normalized wavefunction'
        n_u = mixfunction(u, el)

        'Adding eigenvalue and probability function'
        n_p = mixfunction(p, el)

        plt.plot(x, n_u, 'b')
        plt.plot(x, n_p, 'g')
        plt.plot(x, el, 'r')

    eigen_min = eigen_value/3.0

    if(min(n_u) < y_min):
        y_min = min(n_u)

    if(max(n_p) > y_max):
        y_max = max(n_p)

    'for loop'
    count += 1

'Plotting graph'
plt.figure(1)
plt.xlabel('Position x')
plt.ylabel('Magnitude')
plt.title(title)
plt.legend(['$V(x)$', '$\phi (x)$', '$Probability(x)$', '$E_n$'], loc=4)
plt.axis([potential_initial_position, potential_final_position,
          y_min - 0.1, y_max + 0.1])

plt.show()
