# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Trapping Frequency Variation during Transport with Additional Errors

# <codecell>

# import several packages:
from pylab import *
from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.qubit import Qubit
from sympy.abc import theta, epsilon, phi
from sympy import collect
import panna as p
import numpy as np
from math import pi
from sympy import exp

# <markdowncell>

# $$ |l_{ges}\rangle =  |l_1 \rangle + e^{i\theta} |l_2\rangle + |l_{back}\rangle$$
# "Unitary" Error of the Quarterwave plate ---> In reality $|l_1\rangle$ will not be a perfect $|\sigma^+\rangle$ lattice, but will contain a $|\sigma^-\rangle$ part, the same holds for $|l_2\rangle$. 
# The backrunning beam has double intensity which leads to a factor $\sqrt{2}$ in the electrical field
# $$ |l_1 \rangle = (\cos\epsilon|\sigma^+\rangle - \sin\epsilon |\sigma^-\rangle)e^{ikx} \qquad |l_2\rangle = (\sin\epsilon|\sigma^+\rangle + \cos\epsilon|\sigma^-\rangle)e^{ikx} \qquad |l_{back}\rangle = \sqrt{2}\frac{1}{\sqrt{2}} (|\sigma^+\rangle + |\sigma^-\rangle) e^{-ikx}$$
# For lattice crosstalk occuring due to unorthogonality the resulting states can be written as 
# $$|l_1\rangle = |\sigma^+\rangle e^{ikx} \qquad |l_2\rangle = (\sin\varphi |\sigma^+\rangle + \cos\varphi |\sigma^-\rangle)e^{ikx}$$
# 
# With the states chosen right now, even for perfect $|l_1\rangle,\ |l_2\rangle$ one state will see a contamination from the other lattice:
# $$ U_{|+\rangle} = I_{|\sigma^+\rangle} \qquad U_{|-\rangle} = \frac{7}{8} I_{|\sigma^-\rangle} + \frac{1}{8} I_{|\sigma^+\rangle}$$ 
# sensitivity of the two states to sigma+, sigma-:
# 
# $$ \hat s_0 = 0\cdot |\sigma^+\rangle \langle \sigma^+ | + 1\cdot |\sigma^-\rangle \langle \sigma^-|
# \quad
# \hat s_1 = \sqrt{\frac{7}{8}} \cdot |\sigma^+\rangle \langle \sigma^+ | + \sqrt{\frac{1}{8}}\cdot |\sigma^-\rangle \langle \sigma^-|$$

# <codecell>

sigma_plus = Qubit(0)
sigma_minus = Qubit(1)

s_0 = 0*sigma_plus*Dagger(sigma_plus) + 1*sigma_minus*Dagger(sigma_minus)
s_1 = sqrt(7/8)*sigma_plus*Dagger(sigma_plus)+sqrt(1/8)*sigma_minus*Dagger(sigma_minus) 
s_0, s_1

# <codecell>

# Declaration of symbol names
l1, l2, e1, e2 = Symbol(r'|{l_1}\rangle'), Symbol(r'|{l_2}\rangle'), Symbol('\epsilon_1'), Symbol('\epsilon_2')
k = Symbol('k', real = True)
x = Symbol('x', real = True)
phi = Symbol('phi', real = True)
theta = Symbol('theta', real = True)
l_back = Symbol(r'|l_{back}\rangle')
l_ges = Symbol(r'|l_{ges}\rangle')

alpha, beta, gamma, delta = Symbol('alpha', real = True), Symbol('beta', real = True),Symbol('gamma', real = True),Symbol('delta', real = True)

# <markdowncell>

# For the general case $|l_1\rangle$, $|l_2\rangle$ are written as superposition of $|\sigma^+\rangle$ and $|\sigma^-\rangle$ with coefficients $\alpha$, $\beta$, $\gamma$, $\delta$:

# <codecell>

l1 = (alpha*sigma_plus+beta*sigma_minus)*exp(I*k*x)
l2 = (gamma*sigma_plus + delta*sigma_minus)*exp(I*k*x)
l_back = (sigma_plus + sigma_minus)*exp(-I*k*x)
l_ges = exp(-I*theta)*l1 + exp(I*theta)*l2 + l_back
l1, l2, l_back

# <markdowncell>

# 
# potential seen by state $|0\rangle$ and state $|1\rangle$: 
# $$ U_{|0\rangle} = \langle l_{ges} | \hat s_0^2 |l_{ges} \rangle \quad U_{|1\rangle} = \langle l_{ges} | \hat s_1^2 |l_{ges} \rangle$$

# <codecell>

E0_u = qapply(s_0*l_ges)
I0_u = qapply(Dagger(E0_u)*E0_u)

E1_u = qapply(s_1*l_ges)
I1_u = qapply(Dagger(E1_u)*E1_u)

I0 = powsimp(I0_u).expand(complex = True)
I1 = powsimp(trigsimp(I1_u)).expand(complex = True)

I0

# <markdowncell>

# The intensity distributions $U_{|0\rangle}$ and $U_{|1\rangle}$ contain spatially varying terms $A\cdot \cos(2kx + \theta)+B\cdot \cos(2kx-\theta)$. Using clever trigonometric identities like this one http://en.wikibooks.org/wiki/Trigonometry/Simplifying_a_sin%28x%29_%2B_b_cos%28x%29 one can simplify the spatially varying terms and obtain the effective amplitude:

# <markdowncell>

# <img src="files/Calc2.png" width=600 />

# <markdowncell>

# <img src="files/Calc.png" width=600 />

# <markdowncell>

# <img src="files/photodiode_Fig_1.png" width=400 />

# <codecell>

# the effective amplitude from the equation above can be calculated with the coefficients A and B as:
def Effective_Amplitude(A, B):
    return 0.5*sqrt(simplify((A+B)**2*(cos(theta))**2+(A-B)**2*(sin(theta))**2))

# <codecell>

# collect terms for cos(2*k*x-theta) and cos(2*k*x+theta) to obtain the coefficients A and B
A0 = collect(I0, {cos(2*k*x-theta), cos(2*k*x+theta)}, evaluate = False)
# calculate effective amplitude
A0eff = Effective_Amplitude(A0[cos(2*k*x-theta)], A0[cos(2*k*x+theta)])
A0[cos(2*k*x-theta)],A0[cos(2*k*x+theta)], A0eff

# <codecell>

# the same for state |1>
A1 = collect(I1, {cos(2*k*x-theta), cos(2*k*x+theta)}, evaluate = False)
A1eff = Effective_Amplitude(A1[cos(2*k*x-theta)], A1[cos(2*k*x+theta)])
A1[cos(2*k*x-theta)], A1[cos(2*k*x+theta)],A1eff

# <headingcell level=3>

# Without Errors

# <markdowncell>

# Compare this to page 135, eq. 4.31a and 4.31b in M. Karski's thesis:

# <codecell>

P0 = A0eff.subs({alpha:1, delta:1, beta:0, gamma:0})
P1 = simplify(A1eff.subs({alpha:1, delta:1, beta:0, gamma:0}))
P0, P1

# <headingcell level=3>

# Considering Unorthogonality

# <codecell>

U0 = A0eff.subs({alpha:cos(phi), delta:1, beta:sin(phi), gamma:0})
U1 = A1eff.subs({alpha:cos(phi), delta:1, beta:sin(phi), gamma:0})
U0, U1

# <headingcell level=3>

# Considering Quarterwave plate Errors

# <codecell>

L0 = A0eff.subs({alpha:cos(epsilon), delta:cos(epsilon), beta:-sin(epsilon), gamma:sin(epsilon)})
L1 = A1eff.subs({alpha:cos(epsilon), delta:cos(epsilon), beta:-sin(epsilon), gamma:sin(epsilon)})
trigsimp(L0), trigsimp(L1)

# <headingcell level=3>

# Finding the 0.5 % Region

# <codecell>

def Make_Array(xvals, expr, sub):
    res = []
    for xval in xvals:
        result = real(expr.subs({sub:xval}))
        res.append(result)
    return res

# <headingcell level=4>

# Without Additional Errors

# <codecell>

figure()
p.colorcycle(3)
p.grey_title('Variation of the Trapping Frequency with Errors')
thetas = linspace(0,pi,100)

P = Make_Array(thetas, P0, theta)
plot(thetas, array(P)**0.5, lw = 2, label = r'Trapping Frequency for state $|0\rangle$')

P = Make_Array(thetas, P1, theta)
plot(thetas, array(P)**0.5, '--', lw = 2, label = r'Trapping Frequency for state $|1\rangle$')

#ylim(0.86, 1.005)
xticks([0,pi/2,pi], ['0', r'$\pi/2$', r'$\pi$'])
xlabel(r'Relative Angle $\theta$')
ylabel(r'Trapping Frequency in Multiples of $\omega_{trap}$')

leg = legend()
leg.get_frame().set_alpha(0.4)
#grid()
savefig("Influence_of_Quarterwaveplate_Error_Trapping_Frequency_During_Transport.png")
show()

# <headingcell level=3>

# With Additional Errors

# <codecell>

figure()
p.colorcycle(2)
p.grey_title('Variation of the Trapping Frequency with Errors')
thetas = linspace(0,pi,100)

P = Make_Array(thetas, P1, theta)
phi0 = deg2rad(0.5)
U = Make_Array(thetas, U1.subs({phi:phi0}), theta)
epsilon0 = deg2rad(0.4)
L = Make_Array(thetas, L1.subs({epsilon:epsilon0}), theta)

fill_between(thetas, array(P)**0.5+0.005*array(P)**0.5, array(P)**0.5-0.005*array(P)**0.5, color = 'blue', alpha = 0.4)
plot(thetas, array(L)**0.5, lw = 2.5)
plot(thetas, array(U)**0.5, '--', lw = 2.5)

P = Make_Array(thetas, P0, theta)
phi0 = deg2rad(0.5)
U = Make_Array(thetas, U0.subs({phi:phi0}), theta)
epsilon0 = deg2rad(0.4)
L = Make_Array(thetas, L0.subs({epsilon:epsilon0}), theta)

fill_between(thetas, array(P)**0.5+0.005*array(P)**0.5, array(P)**0.5-0.005*array(P)**0.5, color = 'blue', alpha = 0.4, label = '0.5 % Region for no Error')
plot(thetas, array(L)**0.5, lw = 2.5, label = r'$\epsilon = 0.4^\circ$')
plot(thetas, array(U)**0.5, '--', lw = 2.5, label = r'$\phi =  0.5^\circ$')

ylim(0.86, 1.005)
xticks([0,pi/2,pi], ['0', r'$\pi/2$', r'$\pi$'])
xlabel(r'Relative Angle $\theta$')
ylabel(r'Trapping Frequency in Multiples of $\omega_{trap}$')

leg = legend(loc = 4)
leg.get_frame().set_alpha(0.4)
grid()
show()

# <codecell>


