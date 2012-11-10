# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Phase Noise Measurement via PFD

# <markdowncell>

# If two input signals to the Phase Frequency Discriminator are frequency locked, the error voltage provided by the PFD is proportional to their phase difference. 
# This can be used to obtain the phase noise spectral density
# 
# $$ S_\phi(f) = | \mathcal F \lbrace \phi (t) \rbrace |^2 $$
# 
# by measurement of the power spectral density of the output signal with a spectrum analyzer
# 
# $$ S_V(f) = | \mathcal F \lbrace V(t) \rbrace |^2 . $$
# 
# One advantage of this method compared to the beat signal method is that no carrier is present since for zero phase difference the voltage drops to zero as well. This avoids large requirements on the dynamic range of the spectrum analyzer. Furthermore, the output signal of the PFD shows only phase fluctuations between the signals and has no sensitivity to amplitude fluctuations. 
# 
# On the other hand the measurement takes now place in the region around 0 Hz of the spectrum analyzer, where the device has a large intrinsic $1/f$ noise peak.
# 
# Hence in this method a low noise preamplifier $\triangleright$ is used. To obtain a meaningful result the noise floor of the spectrum analyzer and preamplifier and the amplification factor $\triangleright$ have to be taken into account in the analysis.
# 
# The amplification constant of the low noise preamplifier can be measured by e.g. generating a low frequency (~200 kHz) sinusodial signal and observing one time the amplified version and one time the unamplified version on the spectrum analyzer. The amplification can then be read of in dB.
# 
# Measurement of the system noise floor is done by termination of the amplifier input with the input impedance 50 $\Omega$ of the amplifier. The noise spectrum is recorded, so that it can be substracted from the measured phase noise spectrum. 
# 
# $$ S_V(f) = | \mathcal F \lbrace \phi(t) \cdot K_\phi \cdot \triangleright \rbrace |^2 = K_\phi^2 \cdot \triangleright^2 \cdot |\mathcal F \lbrace \phi(t) \rbrace|^2 .$$
# 
# To obtain the power spectral density with a high resolution and normalized to a 1 Hz bandwidth, several spectra with different spans are recorded and combined as shown in section \ref{SECTION}.
# 
# Finally the rms phase noise $\langle \phi \rangle$ is obtained by considering that it is connected to the spectral distribution by integration \footnote{Actually the filter used by the spectrum analyzer is of a Gaussian form with a width given by the resolution bandwidth. The peak at 0 Hz will therefore "mask" the phase noise in this region, i.e. the spectral distribution will always be the one of the filter shape. Taking this into account, the integration does not involve the lower limit at 0 Hz, but starts from 1.2 Hz.}
# 
# $$ \langle \phi \rangle =\sqrt{ \int_{-\infty}^\infty S_\phi(f) df  }=\sqrt{ 2 \cdot \int_0^\infty  S_\phi(f) df }.$$
# 
# Since the PFD measures the phase difference between the two input channels, the rms deviation $\langle \phi \rangle$ is given by contributions of local oscillator and reference:
# 
# $$ \langle \phi \rangle = \sqrt{ \langle \phi_{LO}\rangle^2 + \langle \phi_{RF} \rangle^2 }.$$

# <codecell>

### The files containing the data for the joined function from the system floor and phase noise measurement
# need to be placed in the same folder as the .py or .ipynb file

# import package for plotting:
from pylab import * 
# import package used for read out of files contained in folder:
import panna as p 
# show indices by which all *.csv files in the folder can be accessed
p.show_n()

# <codecell>

### Insert values for PFD- and amplication-constant and indices for data files###

# Phase Detector Constant of 22 mV / degree in Volt/rad:
K_phi = 180/pi*273e-3 
# amplification in dB
amplifier = 40.0 
# index of the data file that contains the system floor (read off from output from last cell)
ifloor = 1
# index of the data file that contains the phase noise measurement (read off from output from last cell)
iphase = 0

# <codecell>

### Calculation of Phase Noise Spectrum after substracting the floor, amplifier, K_phi
y_res	= 10*log10(abs(10**(array(p.y[iphase])/10.) - 10**(array(p.y[ifloor])/10.))) - 20*log10(K_phi)-amplifier
y_res_lin = 10**(y_res/10.)
x_res	= p.x[iphase]

# <codecell>

n = p.number_of_files()

# integration limits: see footnote in text
xmin = 1.2 # Hz
xmax = 1e6 # Hz

x, y, y_lin, start, stop, phi = [0]*n,[0]*n,[0]*n,[0]*n,[0]*n,[0]*n

# <codecell>

########## Step 1: plot all found input data ############
figure(1)
p.colorcycle(n+1)
for i in arange(0,n): 
	plot(p.x[i], p.y[i], label = str(i))

plot(x_res, y_res, label = r"$S_\phi\ [dB rad^2/Hz]$")

legend()
xlabel("offset frequency [Hz]") 
ylabel("amplitude [dBm/Hz]")
p.grey_title("comparison_step_1_plot_input_data")
savefig("comparison_step_1_plot_input_data.png")
show()

# <codecell>

######### Step 2: find integration limits ############
for i in arange(0,n): 
	######## on linear range #############
	y_lin[i] = [pow(10, element/10.0) for element in p.y[i]] 
	######## find integration limits ############
	start[i] 	= p.find_nearest(p.x[i], xmin)
	stop[i] 	= p.find_nearest(p.x[i], xmax)	

y_res_lin = [pow(10, element/10.0) for element in y_res] 

######## find integration limits ############
start[i] 	= p.find_nearest(x_res, xmin)
stop[i] 	= p.find_nearest(x_res, xmax)

center	= 0

# <codecell>

####### Step 3: Find the phase noise by integration ##########
for i in arange(0,n): 
	### 
	phi[i]	= sqrt(2*abs(np.trapz(y_lin[i][start[i]:stop[i]], p.x[i][start[i]:stop[i]])))
	print str(i) 
	print "power spectral density ---> phase noise = ",round(phi[i],5),"rad (",round(phi[i]/pi*180.0,3),"deg)"

phi_res = sqrt(2*abs(np.trapz(y_res_lin[start[i]:stop[i]], x_res[start[i]:stop[i]])))
print "S_phi (Resulting function):" 
print "power spectral density ---> phase noise = ",round(phi_res,5),"rad (",round(phi_res/pi*180.0,3),"deg)"

### put output into latex form ###

FILE_results = open('results.tex', "w")

results = []
i = 0
for i in arange(0,n):
	results.append(str(i))
	results.append("	power spectral density ---> phase noise = "+str(round(phi[i],5))+"rad ("+str(round(phi[i]/pi*180.0,3))+"deg)")
	results.append(" ")

results.append("S_phi_res")
results.append("	power spectral density ---> phase noise = "+str(round(phi_res,5))+"rad ("+str(round(phi_res/pi*180.0,3))+"deg)")
results.append(" ")

for item in results:
	print >> FILE_results, item

# <codecell>

### Step 4: plot in logscale ###
figure(2)
p.colorcycle(n+1)
for i in arange(0,n): 
	plot(p.x[i], p.y[i], label = str(i))

plot(x_res, y_res, label = r"$S_\phi\ [dB rad^2/Hz]$")

xscale("log")
xlim(0.6, 2.5e6)

grid(True, which = 'both')
legend()
xlabel("Offset Frequency [Hz]") 
ylabel(r"Phase Noise [dB rad$^2$/Hz]")
p.grey_title("Phase Noise Spectrum")
savefig("comparison_step_4_plot_input_data.png")
show()

#### Step 5: Output the resulting S_phi_res function into a file
FILE_my_joined = open('PFD_ResultSPHI.csv', "w")
my_try = []
for i in range(0,len(x_res)):
	my_try.append(str(x_res[i])+","+str(y_res[i]))

for item in my_try:
	print >> FILE_my_joined, item

# <codecell>


