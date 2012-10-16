# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Measurement of the Quarterwave Plate Retardance Error

# <codecell>

import panna as pac
from numpy import * 
from matplotlib import *
from pylab import *

# <headingcell level=3>

# Jones Matrices for Quarterwave Plate with Error & Polarizer

# <codecell>

# Rotation matrix
R = lambda theta: matrix([ [ cos(theta) , -sin(theta) ],
			  [ sin(theta) , cos(theta) ]])

# quarter wave plate with error in retardance
Q_error_p = lambda error: matrix([ [1. , 0.], 
			  [0. , exp(-1j*(pi/2.+error))]])
    
# Linear polarizer with axis of transmission horizontal
P = matrix([ [1. , 0.],
	    [0. , 0.]])

# <headingcell level=3>

# Incoming Electrical Field:

# <codecell>

# incoming field perfectly prepolarized by GP
omega = 80.e6
del_omega = 40.e3
t = linspace(0, 400e-6, 500)
E1 = lambda time_el: exp(1j*(omega+del_omega)*time_el)*matrix([ 	[1.],
	     						[0.]])

# <headingcell level=3>

# Intensities on Photodiode with and without Quarterwave Plate:

# <codecell>

# without quarter wave plate
def signal(t, angle): 
	I = []
	for el in t:
		E = R(angle)*P*R(-angle)*E1(el)
		I.append(E.conj().T*E)
	return array(I)[:,0]

# with imperfect quarter wave plate 
def Q_error_plus(t, angle, error): 
	I = []
	for el in t:
		E = R(angle)*P*R(-angle)*R(pi/4.)*Q_error_p(error)*R(-pi/4.)*E1(el)
		I.append(E.conj().T*E)
	return array(I)[:,0]

# <headingcell level=3>

# Applying above Functions and Plotting them:

# <codecell>

thetas_deg = linspace(0, 1200, 300)
thetas_rad = thetas_deg*2*pi/360.

V_without	= []
V_perfect 	= []
V_plus 		= []

for theta in thetas_rad:
	V_without.append(real(signal([0.0], theta)[0]))
	V_perfect.append(real(Q_error_plus([0.0], theta, 0)[0]))
	V_plus.append(real(Q_error_plus([0.0], theta, 2*pi/40.)[0]))

def angle(time):
	scale = 180./0.15 # scaling factor resulting from motor speed/ rotation speed of the polarizer
	return (array(time)-time[0])*scale

def colorcycle():
	colormap = plt.cm.gist_ncar
	num_plots = 5
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])

# <codecell>

figure()

one = subplot(311)
one.set_title("expectation", style='italic', color = 'white',
        bbox={'facecolor':'grey', 'alpha':0.5})
colorcycle()
color1 = cm.gist_ncar(0.75)

one.set_xticks([])
one.set_xlim(0,1200)
one.set_ylim(0,1)
one.set_yticks([0,1])

plot(thetas_deg, V_without, label = r"without $\lambda/4$ plate")
plot(thetas_deg, V_perfect, label = r"$\lambda/4$ perfect")
plot(thetas_deg, V_plus, label = r"$\lambda/4 \pm \lambda/40$")

legend()
two = subplot(312)
two.set_title("experimental data", style='italic', color = 'white',
        bbox={'facecolor':'grey', 'alpha':0.5})
colorcycle()
two.set_xticks([])
two.set_xlim(0,1200)
two.set_ylim(0,1)
two.set_yticks([0,1])

plot(angle(pac.x[0]), array(pac.y[0])/max(pac.y[0]), label = r"without $\lambda/4$ plate")
plot(angle(pac.x[1])-40., array(pac.y[1])/max(pac.y[0])+1/2.,'--',lw =2, label = r"with $\lambda/4$ plate", color = color1)

legend()

#### 
three = subplot(313)
three.set_title("zoom into the wiggling", style='italic', color = 'white',
        bbox={'facecolor':'grey', 'alpha':0.5})
colorcycle()

wiggling = array(pac.y[1])+0.5
three.set_xlim(0,1200)
three.set_ylim(0.4,0.6)
three.set_xlabel(r"angle ($^\circ$)")

V_t, V_40, V_100, V_320, V_200    = [], [], [], [], []

for theta in thetas_rad:
    V_t.append(real(Q_error_plus([0.0], theta, 2*pi/12.)[0]))
    V_40.append(real(Q_error_plus([0.0], theta, 2*pi/40.)[0]))
    V_100.append(real(Q_error_plus([0.0], theta, 2*pi/150.)[0]))
    V_320.append(real(Q_error_plus([0.0], theta, 2*pi/320.)[0]))
    V_200.append(real(Q_error_plus([0.0], theta, 2*pi/200.)[0]))

plot(thetas_deg, V_t, label = r"$\lambda/4 \pm \lambda/12$")
plot(thetas_deg, V_40, label = r"$\lambda/4 \pm \lambda/40$")
plot(thetas_deg, V_100, label = r"$\lambda/4 \pm \lambda/150$")
plot(thetas_deg, V_320, label = r"$\lambda/4 \pm \lambda/320$")

plot(angle(pac.x[1])-80., wiggling,'--', lw =2,   color = color1, label = r"data with $\lambda/4$ plate")

legend()

tight_layout()
savefig('lambda_4_plate_error.png')
show()

print "Wiggling Amplitude =", max(wiggling)-min(wiggling)
print "L150 Amplitude =", max(V_100)-min(V_100)
print "L320 Amplitude =", max(V_320)-min(V_320)
print "L200 Amplitude =", max(V_200)-min(V_200)

# <headingcell level=2>

# ##############################################################################################

# <headingcell level=1>

# Percentage of circular polarization with retardance error

# <codecell>

from pylab import *
#from numpy import *
from sympy import *
from sympy.physics.quantum import Dagger
from sympy.abc import theta, epsilon
import panna as p
import numpy as np
from math import pi
import sympy 

# <markdowncell>

# Polarization basis vectors:

# <codecell>

lin_x = Matrix([[1],[0]])
lin_y = Matrix([[0],[1]])
sigma_plus  = 1/sqrt(2) * Matrix([[1],[1j]])
sigma_minus = 1/sqrt(2) * Matrix([[1],[-1j]])
k = 2*pi/(865.9e-9)

# <codecell>

# Rotation matrix:
theta = Symbol('theta')
R = lambda theta: Matrix([[cos(theta), sin(theta)],[-sin(theta), cos(theta)]])
R(theta)

# <headingcell level=4>

# Jones Matrix for a rotated waveplate with retardance error

# <codecell>

# waveplate with retardance epsilon and error rotated to some angle:
class waveplate:
    def __init__(self, epsilon, error, angle):
        self.epsilon = epsilon
        self.error = error
        self.angle = angle
    def M(self):
        # http://en.wikipedia.org/wiki/Jones_calculus 
        # retardance epsilon = pi for lambda/2; pi/2 for lambda/4
        return R(-self.angle)*Matrix([[1,0],[0,sympy.exp(I*(self.epsilon+self.error))]])*R(self.angle)
    def update_error(self, new_error):
        self.error = new_error
    def update_angle(self, new_angle):
        self.angle = new_angle

# del_x is the waveplate error in multiples of lambda, returned is the error in rad
def err(del_x): 
    return del_x*2*pi

def Norm(some):
    return abs(complex(some[0]))**2

# <markdowncell>

# Often the retardance error is expressed in multiples of $\lambda$. The retardance error in rad is calculated by using (--> err(del_x))
# $$ \frac{\Delta x}{\lambda} = \frac{\Delta\phi}{2\pi}.$$

# <markdowncell>

# Test of the ansatz for the quarterwave plate: 
# 
# For a quarterwave plate rotated to a 45° (pi/4) angle without retardance error one should obtain circular polarization from linear incoming polarization. 

# <codecell>

wp4 = waveplate(pi/2, 0, pi/4)
wp4.M(), wp4.M()*lin_x, Norm(Dagger(sigma_plus)*wp4.M()*lin_x), Norm(Dagger(sigma_minus)*wp4.M()*lin_x)

# <markdowncell>

# The output numbers below show the percentage of circular polarization, resulting from the projection on $|\sigma^+\rangle$ and $|\sigma^-\rangle$ for linearly polarized input states:
# $$ |\langle \sigma^- | \Lambda_4(\epsilon, \pi/4) | x\rangle|^2 $$ 

# <codecell>

wp4 = waveplate(pi/2, 0, pi/4)
wp4.M()

Norm(Dagger(sigma_minus)*wp4.M()*lin_x),Norm(Dagger(sigma_plus)*wp4.M()*lin_x), Norm(Dagger(sigma_minus)*wp4.M()*lin_y), Norm(Dagger(sigma_plus)*wp4.M()*lin_y)

# <codecell>

wp4.update_error(err(1/450))

Norm(Dagger(sigma_minus)*wp4.M()*lin_x),Norm(Dagger(sigma_plus)*wp4.M()*lin_x), Norm(Dagger(sigma_minus)*wp4.M()*lin_y), Norm(Dagger(sigma_plus)*wp4.M()*lin_y)

# <codecell>

wp4.update_error(err(1/200))

Norm(Dagger(sigma_minus)*wp4.M()*lin_x),Norm(Dagger(sigma_plus)*wp4.M()*lin_x), Norm(Dagger(sigma_minus)*wp4.M()*lin_y), Norm(Dagger(sigma_plus)*wp4.M()*lin_y)

# <codecell>

wp4.update_error(err(1/12))

Norm(Dagger(sigma_minus)*wp4.M()*lin_x),Norm(Dagger(sigma_plus)*wp4.M()*lin_x), Norm(Dagger(sigma_minus)*wp4.M()*lin_y), Norm(Dagger(sigma_plus)*wp4.M()*lin_y)

# <codecell>

wp4.update_error(err(1/150))

Norm(Dagger(sigma_minus)*wp4.M()*lin_x),Norm(Dagger(sigma_plus)*wp4.M()*lin_x), Norm(Dagger(sigma_minus)*wp4.M()*lin_y), Norm(Dagger(sigma_plus)*wp4.M()*lin_y)

# <codecell>


