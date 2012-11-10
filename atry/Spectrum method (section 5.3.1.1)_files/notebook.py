# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Phase Noise via Beat Signal

# <markdowncell>

# 
# The spectrum of an ideal waveform $\propto\sin(\omega_0 t)$ is a $\delta$-peak at the carrier frequency $\omega_0$ which contains all the signal power. Now, if one supposes a phase modulation $\phi (t)$ the resulting signal is
# $$ \mathcal S_{ideal}(t) = \sin(\omega_0 t) \Rightarrow S_{ideal}(f)=\delta(\omega_0) \qquad \mathcal S_{mod}(t) = \sin(\omega_0 t + \underbrace{ m \sin(2\pi f_m t) }_{\phi(t)}.$$
# 
# For this the mean square phase deviation can be calculated as 
#     $$\langle \phi \rangle^2 = \lim_{T\rightarrow \infty} \frac{1}{T} \int_{-T/2}^{T/2} |\phi(t)|^2 dt = \frac{m^2}{2} .$$
# 
# The modulated spectrum $\mathcal S_{mod}(f)$ is described by Bessel functions, which appear as sidebands to the carrier at multiples of the modulation frequency $f_m$.
# If the modulation is small enough the higher order Bessel coefficients become negligibel and the modulated signal can be written as
# $$ e^{i( \omega t + m\sin(\omega_m t))} = e^{i\omega t} \sum_{n=-\infty}^\infty J_n(m) e^{in \omega_m t} \stackrel{m < 0.4}{\approx}  e^{i\omega t} \left[ \frac{m}{2} e^{-i\omega_m t} + 1 + \frac{m}{2}e^{i\omega_m t} \right]. $$
# 
# [FIGURE BELOW ABOUT THIS LIMIT]
# 
# In this case two sidebands with power $m^2/4$ in each will appear. This power will be missing in the carrier. The amplitudes of the sidebands and the carrier are connected via the modulation index by 
# 
# $$ \frac{P_{single\ sideband}}{P_{carrier}} = \frac{m^2}{4}= \frac{\langle \phi \rangle^2}{2}$$
# 
# (see e.g. section \ref{SECTION ABOUT BESSEL FUNC}).
# 
# A continous noise spectrum can be "build up" by small modulations at different frequencies. The spectrum on the spectrum analyzer is measured in dBm at a certain resolution bandwidth RBW in Hz. The measured power is normalized to a 1 Hz by substracting $10 \cdot log_{10}(RBW)$. \\
# 
# From this one can find a connection between the mean square phase deviation and integration over the normalized signal spectrum seen on the spectrum analyzer:
# 
# $$ \langle \phi \rangle^2 = \frac{ \int_{both\ sidebands} P(\nu)[dBm/Hz] d\nu }{ P_{carrier} } = \frac{ 2 \int_{one\ sideband} P(\nu) d\nu }{P_{carrier}} = 2 \int_{one\ sideband} P(\nu)[dBc/Hz] d\nu.$$
# 
# The unit dBc means that this spectrum has been normalized to the carrier amplitude.
# 
# Alternatively this can be written in a "power fraction in carrier"-form used e.g. in PREVEDELLI
#     $$ \eta = 1 - \langle \phi \rangle^2 = \frac{ P_{carrier} }{ \int_{-\infty}^\infty P(\nu) [dBm/Hz] d\nu} $$
# 
# Advantages of this method include that no quadrature has to be kept, which allows also for the measurement of relatively large phase noise or signals that drift in phase during the measurement time. The limit in which this approximation is still valid ($m<0.4$) means that the suppression of phase noise signal compared to the carrier has to be larger than -14 dB \footnote{This is $10 \log (0.4^2/4)$.} for all frequencies. Still the frequency stability has to be given, such that the carrier frequency stays at the same value. \\
# 
# On the other hand this method is sensitive to amplitude fluctuations. To obtain the rms-deviation most spectrum analyzers offer an "I-O measurement" \footnote{Available for the used Agilent N9010A for several thousand Euros, but still limited by the same dynamic range.} option to differ between the noise which is in-phase with the carrier (amplitude) and in-quadrature with the carrier (phase-noise). 
# 
# In the case of the used external-cavity diode lasers the amplitude noise is negligible compared to the noise added by elements in the optical path. For WHICH AND WHY lasers this is not the case. \\
# If the phase noise is decreased by the phase locked loop, the amplitude noise can no longer be neglected. \\
# 
# Another limit to this measurement method is the free dynamic range of the spectrum analyzer. The free dynamic range at a certain offset frequency gives the smallest signal that can be displayed together with a carrier of a certain strength. In each measurement it has to be verified that the spectrum signal lies several dB above the spectrum analyzers dynamic range. \\

# <headingcell level=3>

# Modulation Index for which higher orders are negligible

# <codecell>

# The joined and RBW normalized function is contained in a file in the folder and was calculated by the "Joined-Function Code
import panna as p # loads data contained in folder
import numpy as np		# numpy for array handling, integration
from numpy import array, arange
from math import sqrt, log10 ,pi, log
from pylab import *
from scipy.special import jn # jn(n,x) returns the Bessel function of integer order n at x

# <codecell>


mod_ind = linspace(0,10,100)
figure(1)
subplot(121)
for i in arange(0,3):
    plot(mod_ind, jn(i, mod_ind), label=r"$J_{"+str(i)+"}(m)$")
legend()
xlabel(r"$m$")
ylabel(r"Bessel Function Value $J_n (m)$")

mod_ind = linspace(0,1,100)
subplot(122)
for i in arange(0,3):
    plot(mod_ind, jn(i, mod_ind), label=r"$J_{"+str(i)+"}(m)$")
plot(mod_ind, mod_ind/2, '--', label = r"$A(m) = m/2$")

xlabel(r"$m$")
ylabel(r"Bessel Function Value $J_n (m)$")
legend()
tight_layout()
savefig("Phase_Noise_via_Beat_Signal_Bessel_functions.png")
show()

# <markdowncell>

# If the modulation is small the higher order Bessel coefficients become negligibel and only two sidebands with power $m^2/4$ in each will appear. This power will be missing in the carrier. The amplitudes of the sidebands and the carrier are connected via the modulation index by 
#     $$ \frac{P_{single\ sideband}}{P_{carrier}} = \frac{m^2}{4}= \frac{\langle \phi \rangle^2}{2}.$$

# <codecell>

# Integration limits
xmin = 1.57 # Hz, otherwise the shape of the resolution bandwidth is interpreted as phase noise
xmax = 1e6 # Hz

# <codecell>

data = p.load_all()
len_files = p.number_of_files()
print len_files

x, y = [0]*len_files, [0]*len_files
for i in arange(0,len_files):
    x[i], y[i] = list(data[i][:,0]), list(data[i][:,1])

# <headingcell level=2>

# Plot found Files

# <codecell>

figure(1)
p.colorcycle(len_files)
for i in arange(0, len_files): 
	plot(x[i], y[i], label = str(i))
	i += 1

legend()
xlabel("offset frequency [Hz]") 
ylabel("amplitude [dBm/Hz]")
p.grey_title("comparison_step_1_plot_input_data")
savefig("comparison_step_1_plot_input_data.png")
show()

# <headingcell level=2>

# Shift to 0 Hz Center

# <codecell>

def find_center(x_values, y_values):
	y_values = list(y_values)
	maximum_index 	= y_values.index(max(y_values))	
	centerfreq	= x_values[maximum_index]
	return centerfreq, maximum_index

for i in arange(0, len_files):
    x[i] = list(array(x[i])-find_center(x[i],y[i])[0])

# <headingcell level=2>

# dBc values 

# <codecell>

######### Step 2: dBc values and integration limits ############
y_lin, start, stop = [0]*len_files, [0]*len_files, [0]*len_files

for i in arange(0, len_files): 
	######## dBc and linear values #############
	y[i] = array(y[i]) - max(y[i])
	y_lin[i] = [pow(10, element/10.0) for element in y[i]] 
	######## find integration limits ############
	start[i] 	= x[i].index(min(x[i], key= lambda d: abs(xmin-d)))
	stop[i] 	= x[i].index(min(x[i], key= lambda d: abs(xmax-d)))	

center	= 0

# <headingcell level=2>

# Result / Integration

# <markdowncell>

# Equation 1:    
# $$ \langle \phi \rangle^2 = \frac{ \int_{both\ sidebands} P(\nu)[dBm/Hz] d\nu }{ P_{carrier} } = \frac{ 2 \int_{one\ sideband} P(\nu) d\nu }{P_{carrier}} = 2 \int_{one\ sideband} P(\nu)[dBc/Hz] d\nu $$
# Equation 2:
# $$ \eta = e^{-\langle \phi\rangle^2} =  1 - \langle \phi \rangle^2 = \frac{ \int P_{carrier} d\nu }{ \int_{-\infty}^\infty P(\nu) [dBm/Hz] d\nu} = \frac{\int_0^{RBW_{min}} P_{carrier} d\nu}{\int_{0}^\infty P(\nu) d\nu }$$

# <codecell>

####### Step 3: Find the phase noise by integration ##########
phi, mu = [0]*len_files, [0]*len_files

for i in arange(0, len_files): 
    ### Equation 1
    phi[i]	= sqrt(2*abs(np.trapz(y_lin[i][start[i]:stop[i]], x[i][start[i]:stop[i]])))
    ### Equation 2
    mu[i]	= np.trapz(y_lin[i][center:start[i]], x = x[i][center:start[i]])/np.trapz(y_lin[i][center:stop[i]], x = x[i][center:stop[i]])
    print "file: ",i 
    print "Eq1: power spectral density ---> phase noise = ",round(phi[i],5),"rad (",round(phi[i]/pi*180.0,3),"deg)"
    print "E12: power in carrier ---> phase noise = ",round(sqrt(-log(mu[i])),5),"rad (",round(sqrt(-log(mu[i]))/pi*180.0,3),"deg)"

# <codecell>

### Step 4: plot in dBc and logscale ###
figure(2)
p.colorcycle(len_files)
for i in arange(0, len_files):
	plot(x[i], y[i], label = r"file "+str(i)+r"; $\langle \phi \rangle $ = "+str(round(phi[i]/pi*180.0,3))+r"$^\circ$")

plot([xmin]*10,linspace(min(y[0]), max(y[0]), 10), '--', color = "grey")
plot([xmax]*10,linspace(min(y[0]), max(y[0]), 10), '--', color = "grey", label = "Integration limits")
xscale("log")
xlim(0.54, 2e6)

legend()
xlabel("offset frequency [Hz]") 
ylabel("amplitude [dBc/Hz]")
p.grey_title("Phase Noise in dBc/Hz")
grid(True, which = 'both')
savefig("comparison_step_4_plot_input_data.png")
show()

# <codecell>

### Step 5: output files in dBc ###
for i in arange(0, len_files):
	print "saving dBc-joined function for index "+str(i)
	my_filename = str(i)+"_dBc.csv"
	FILE_my_filename = open(my_filename, "w")
	i = 0
	my_try = []
	for j in range(0, len(x[i])):
		my_try.append(str(x[i][j])+","+str(y[i][j]))
	for item in my_try:
		print >> FILE_my_filename, item
	i += 1

# <codecell>


