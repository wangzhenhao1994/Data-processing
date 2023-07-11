1D-Schroedinger-Equation
========================

This is a python program to solve 1D Schroedinger Equation (SE) for eigenvlaues and eigenfunctions.

SE has played a very important role for quantum mechanical particles. The behavior of the quantum particles can be understood by solving SE. The stationary state of the particles can be obtained from the SE for the appropriate potential.

This code has been developed as a class work for the student of Computational Physics, Department of Physics, Shahjalal University of Science and Technology, Sylhet - 3114, Bangladesh by Md. Enamul Hoque, Lecturer, Department of Physics.

email: mjonyh-phy@sust.edu or mjonyh@gmail.com
facebook: www.facebook.com/mjonyh

1. How to run?
-------------------------
	$ python schroedinger.py
    or
    $ python2 schroedinger.py
	
2. Used method:
-------------------------
	1. Numerov Algorithm: for boundary value.
	2. Simpson Method: for normalization.
	3. Secant Method: for searching eigenvalue.

3. How to include potential?
------------------------------
	1. In the top of schroedinger.py file, there is a function named 'v(x)' for the potential. If the potential has an analytical form than just put the equation for 'x' in there.
	2. If the potential is in tabular form, it can also be called from 'v(x)'. (Under development)

