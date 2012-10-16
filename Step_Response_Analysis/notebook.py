# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Final Step Response Analysis

# <markdowncell>

# Here the final analysis of the phase locked loop step response contained in the results section of my thesis is obtained. 

# <codecell>

### Definition of some Functions ###
# calculates derivative
def dy_dx(x,y):
    dx = x[1::] - x[:-1] 
    dy = y[1::] - y[:-1]
    timestep = (x[-1]-x[0])/(len(dx))
    return x[:-1]+timestep/2., dy/dx

# Gives Spectrum of Data, Tried out below
def frequency_analysis(time, signal):
    spectrum = np.fft.fft(signal)#*np.hanning(len(signal)))
    dt = abs(time[-1] - time[0])/len(time)
    frate = 1/dt
    freq = np.fft.fftfreq(len(signal))
    freq_in_hertz=abs(freq*frate) # from http://stackoverflow.com/questions/3694918/how-to-extract-frequency-associated-with-fft-values-in-python
    """ Check if the frequencies are in agreement with Nyquist sampling theorem:
    f_min_exp = 1/2*1/(time[-1]-time[0])
    f_max_exp = 1/2*1/dt
    print f_min_exp, f_max_exp, min(freq_in_hertz), max(freq_in_hertz)"""
    print min(freq_in_hertz)
    # Remove small high frequency part
    tol = 0.05*abs(spectrum).max()
    for i in xrange(len(spectrum)-1, 0, -1):
        if abs(spectrum[i]) > tol:
            break
    return freq_in_hertz[:i+1], spectrum[:i+1]

def data_region(x_array, y_array, xmin, xmax):
    x_new = x_array[p.find_nearest(x_array, xmin) : p.find_nearest(x_array, xmax)]
    y_new = y_array[p.find_nearest(x_array, xmin) : p.find_nearest(x_array, xmax)]
    return x_new, y_new

# <codecell>

### Import used Packages ###
import panna as p  # indices shown below are ok
from numpy import *
from pylab import *
from numpy import exp, cos, sqrt, sin

p.show_n()

# <codecell>

### Rename the Data that is obtained from the *.csv-files with "p.column[index from above]" into a
# human readable form
### the data region of the actual step is chosen and a factor is applied, such that 
# all data is normalized between 0 and 1, where 1 is the set value
taom0, raom0 = data_region(p.x[4], p.y[4]*1/0.7, 0.0, 15e-6) 
taom1, raom1 = data_region(p.x[5], p.y[5]*1/0.7, 0.0, 6e-6)
tVCO, rVCO = p.x[6], p.y[6]
tEOM, rEOM = p.x[3], p.y[3]
tDDS, rDDS = p.x[2], p.y[2]

figure()
plot(tVCO, rVCO)
plot(taom0, raom0)
show()

# <headingcell level=2>

# Connection of Step Response and System Transfer Function

# <markdowncell>

# Adjustment of the system parameters for the dynamic behavior observed in the step response is done in a trade off between speed and stability. The system can be brought into a state, where it reacts very fast to a change in the reference, but then it might oscillate around the setpoint, instead of fastly settling on it. 
# 
# The phase noise spectral density gives mainly information about the resulting stability of the system. As suggested a stable or unstable behavior and the BANDWIDTH of the system can be estimated by height and position of the servo-bump.
# 
# An elegant way to observe the dynamic behavior and to extract the system transfer function $H(s)$ at the same time is to record the closed-loop step response: This is the time response of the controlled oscillator phase $\phi_{osc}$ when a phase step is applied to the reference $\phi_{ref}$. Four quantities are used to characterize the closed-loop step response:
# 
# 	* Rise Time: The output rises beyond 90 \% of the desired value for the first time. The rise time should be as short as possible, such that the system responds fast and minimizes the difference to the reference again.
# 	* Overshoot: How much higher is the peak level than the steady state? The error between reference and controlled oscillator will decrease when the oscillator phase approaches the reference level. Since the oscillator phase depends on error signals from longer times (see CONVOLUTION EQ) a phase locked system with fast response and therefore small damping will show an overshoot. 
# 	* Settling Time: The time it takes to converge to the steady state. After the overshoot the system will settle into the steady state error region after the settling time, which should be as small as possible. Another common definition is to define the signal as settled when it stays within the 2 \% or 5\ \% region of the set value.
# 	* Steady State Error: The difference between the steady-state and the desired output for $t\rightarrow \infty$. 
# 
# DO IT YOURSELF SOS:
# [http://www.atp.ruhr-uni-bochum.de/rt1/syscontrol/node57.html#fig:8.2.1] 
# 
# 
# To obtain good parameters for the phase locked loop, a phase step can be applied to the reference and the characteristics above are optimized by changing the available degrees of freedom contained in the system transfer function. 
# 
# From the phase step response the system transfer function can be obtained:
# 
# In the Laplace domain a unit step of the reference signal $\phi_{ref}(s)$ is 
# $$ \phi_{ref}(s) =  \mathcal L \lbrace \phi_{ref}(t) \rbrace \stackrel{unit\ step}{=} \int_0^\infty e^{-st} dt = \frac{1}{s},$$
# 
# which gives a connection of controlled oscillator and system transfer function
# $$ H(s) = \frac{\phi_{osc}(s)}{\phi_{ref}(s)} \quad \Leftrightarrow \quad \phi_{osc}(s) = \frac{1}{s} H(s).$$
# 
# Additionally it can be used that taking the derivative in the time domain leads to multiplication with a factor $s$ in the Laplace domain. 
# $$ H(s) = \mathcal L \lbrace \frac{d}{dt} \phi_{osc}(t) \rbrace (s).$$
# 
# For a stable system with small steady state error the derivative $\frac{d}{dt} \phi_{osc} (t)$ will approach zero for $t \rightarrow \infty$. Therefore the Laplace transform can be replaced by the Fourier transform in this case.

# <headingcell level=2>

# Estimation of the System Transfer Function

# <markdowncell>

# The open loop transfer function $G(s)$, introduced in section \ref{INTRO} is obtained by multiplication of the transfer-functions of the single components. As discussed in section \ref{PID} the proportional gain $K_P$ gives the main contribution to the PID controllers transfer function. The delay caused mainly by the AOM is approximated by \cite[page 312]{Lunze} 
# $$ e^{-sT_{dead}} \approx \frac{1}{T_{dead} +1}.$$
# 
# Therefore the open loop transfer function is approximately given by
# $$ G^{approx.}(s) \approx \frac{ K_\phi K_P K_{VCO} }{s(T_{dead} s + 1)} = \frac{\mathcal K}{s(T_{dead}s+1)},$$
# 
# and the system transfer function can be obtained from the relation 
# 
# $$ H(s) = \frac{G(s)}{1+G(s)} \quad \Rightarrow \quad H^{approx.}(s) = \frac{ \frac{\mathcal K}{T_{dead}} }{s^2 +  2 \cdot \sqrt{ \frac{ \mathcal K }{T_{dead}}} \cdot \frac{1}{2\sqrt{T_{dead} \mathcal K} }s+ \frac{\mathcal K}{T_{dead}}}.$$ 
# 
# The resulting expression for $H^{approx.}(s)$ has the standard form of a system transfer function for a so called second order system \cite[page 176]{Ogata}:
# $$ H^{2nd\ order}(s) = \frac{\omega_n^2}{s^2 + 2 \zeta \omega_n s + \omega_n^2},$$
# 
# i.e. the system and its dynamics are described by two parameters 
# 
# $\zeta$: damping
# and
# $\omega_n$: natural frequency. 
# 
# For a unit step input one obtains 
# 
# \label{eq:Ct}
# $$ C(s) = H(s) \frac{1}{s} \qquad C(t) = \mathcal L^{-1} \lbrace C(s) \rbrace = 1-e^{-\zeta \omega_n t} \left( \cos\omega_d t + \frac{\zeta}{\sqrt{1-\zeta^2}} \sin\omega_d t \right) \qquad \omega_d = \omega_n \sqrt{1- \zeta^2 }.$$

# <headingcell level=2>

# Second Order System Characteristics

# <markdowncell>

# ZETA REGION
# 
# For a stable, damped (non critically or critically overdamped or undamped) system $\zeta \in ]0,1[$ \cite[page 285]{Lunze}.
# 
# MAXIMUM OVERSHOOT
# 
# At the point of the maximum overshoot the function $C(t)$, which describes the step response, reaches its maximum value. One can either calculate the extrem value by taking the derivative of $C(t)$ or argue that for times $\omega_d t= \pi, 3\pi/4$ (and full multiples) the contribution of the $\cos(\omega_d t)$ or $\sin(\omega_d t)$ gives a contribution $-1$ and therefore enhances $C(t)$. Furthermore for values of $\zeta \in ]0,1[$ the $\sin$ prefactor $\zeta/\sqrt{1-\zeta^2}$ is always smaller one, such that the overshoot will occur at $ \omega_d t_{M} = \pi $ (and odd integer multiples).
# 
# The resulting overshoot in percent is 
# [ALIGN]
# $$ M = C(t_{M}) - 1 \stackrel{t_M = \pi/\omega_d}{=} e^{-\zeta \omega_n (\pi/\omega_d)} \left( \cos\pi + \frac{\zeta}{\sqrt{1-\zeta^2}} \sin\pi \right) = e^{-(\zeta/\sqrt{1-\zeta^2} )\pi}$$
# 
# The maximum overshoot depends only on the damping ratio $\zeta$. For a second order system the phase margin $\varphi_{margin}$ introduced in section \ref{INTRO} and the overshoot can be related by 
# $$ \varphi_{margin}[^\circ] + M [\%] \approx 70 $$
# (derivation in appendix).
# Since the phase margin should be between $30^\circ\dots 60^\circ$ the optimum region for the overshoot is between $10\%\dots 40\%$.
# 
# 
# (This result is true also for $\zeta$ outside the assumed region, see \cite[page 284]{Lunze})

# <codecell>

figure("Maximum Overshoot vs Damping")
ax = subplot(111)
zeta = linspace(0,1,100)
plot(zeta, exp(-(zeta/sqrt(1-zeta**2))*pi), lw = 1.5)
#,
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_major_locator(MultipleLocator(0.1))

grid(True, which = 'both')
#yticks(arange(0.0,1.1,0.1), [r'0%', r'10%', '20%', '30%', '40%', r'50%', '60%', '70%', '80%', '90%', r'100%'])
xlabel(r'Damping Ratio $\zeta$')
ylabel(r'Maximum Overshoot $M$ in Percent')
tight_layout()
savefig("Step_Response_Maximum_Overshoot.png")
show()

# <markdowncell>

# SETTLING TIME
# 
# The step response $C(t)$ given in equation \ref{eq:Ct} has an envelope $1-e^{-\zeta \omega_n t}$.
# Since the oscillating part will always stay within the envelope, the settling time can be expressed in terms of the envelope function. To stay within the 2 \% region around the set value the signal needs the settling time (see also \ref[page 183]{Ogata})
# \label{eq:SecondOrderSettle}
# $$ e^{-\zeta\omega_n t_{settle}} \approx 0.02 \quad \Rightarrow \quad t_{settle} \approx \frac{4}{\zeta \omega_n}.$$
# 
# Since both, settling time and overshoot reach better values for larger damping ratio, the optimum region for $\zeta$ lies around $\zeta_{opt} = 0.7$ \cite{}

# <codecell>

figure("Settling Time")

zeta = linspace(0,1,100)
plot(zeta, -1/zeta*log(0.02*sqrt(1-zeta**2)), lw = 1.5)
ylim(0,15)
tight_layout()
show()

# <headingcell level=2>

# Approximation as Second Order System

# <headingcell level=3>

# I) Second Order System Parameters from Step Response

# <codecell>

figure("AOM Response")
ax = subplot(111)
p.colorcycle(3)
p.grey_title("Measured System Step Response")

ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.25e-6))
ax.xaxis.set_major_locator(MultipleLocator(1e-6))

#fill_between(tDDS, 0.99, 1.01, color = 'grey', alpha = 0.75)
amp = linspace(0,1,10)
time = linspace(-1.4e-6, 6.4e-6, 10)
plot([0.260e-6]*len(amp), amp, '--', lw = 1.5, color = 'grey')
plot([0.520e-6]*len(amp), amp, '--', lw = 1.5, color = 'grey')
plot([2e-6]*len(amp), amp, '--', lw = 1.5, color = 'grey')
plot(time , [1.43]*len(time), '--', lw = 1.5, color = 'grey')

grid(True, which = 'both')
plot(tDDS, rDDS, lw = 1.5, color = 'grey', label = r"$C_{Reference} (t)\ (Trigger)$")
plot(taom0, raom0, lw = 1.5, label = r"$C_{Optics}(t)$")

#plot(tVCO, rVCO)
#plot(taom1, raom1)
ylim(0, 1.5)
xlim(-1.4e-6, 6.4e-6)
p.scale_xaxis(ax,1e6)
xlabel(r"Time $\mu$s")
ylabel(r"Amplitude relative to Set Value")
yticks(arange(0.0,1.6,0.2), [r'0%', '20%', '40%', '60%', '80%', r'100%', '120%', '140%'])

legend(loc = 4)
#grid()
tight_layout()
savefig("Final_Step_Response_Optics.png")
show() 

# <headingcell level=3>

# II) Second Order System Parameters expected from Device Characteristics

# <codecell>

def w(Kphi, Kp, KVCO, Tt):
    return sqrt(Kphi*Kp*KVCO/Tt)

def Zeta(Kphi, Kp, KVCO, Tt):
    return 1/(2*sqrt(Tt*Kphi*Kp*KVCO))

Kphi = 1.26 # V/rad
KVCO = 9e6/(2*pi) # rad/(s*V)
Tt = 270e-9
Kp = 8.2
print "omega_n expected =", w(Kphi, Kp, KVCO, Tt),"Zeta expected =", Zeta(Kphi, Kp, KVCO, Tt)

# <headingcell level=3>

# Obtaining the Transfer Function for the AOM

# <markdowncell>

# The unapproximated transfer function from which bandwidth and phase stability can be calculated, is obtained from the step response as introduced in section \ref{SEC}. To investigate the exact influence of the AOM-dead time, the step response of the circuit with AOM is compared to the one of the purely electronic circuit:
# In this case the output signal of the VCO does not drive the AOM, but its output is directly connected to the PFD instead of the beat signal recorded on the photodiode. The result is shown in figure \ref{fig:stepElectronics}. 
# 
# The overshoot is only $ M = 2 \%$ which indicates a damping value of $\zeta \approx 0.78$. This step response shows a short settling time of about 800 ns into the 2 \% set value region. It indicates optimal tuning since the settling time has an absolute minimum at $\zeta = 0.76$ \cite[page 183]{Lunze}. 
# 
# Remarkable is that the optical signal crosses the 90 \% value at the same time with the electronic signal, even though it crosses the 10 \% value 100 ns later. This again shows the low damping of the optical signal compared to the electronic. 

# <headingcell level=4>

# Step Response of Electronics

# <markdowncell>


# <codecell>

figure("AOM Response")
ax = subplot(111)
p.colorcycle(3)
p.grey_title("Measured System Step Response")

ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.25e-6))
ax.xaxis.set_major_locator(MultipleLocator(1e-6))

#fill_between(tDDS, 0.99, 1.01, color = 'grey', alpha = 0.75)
"""
amp = linspace(0,1,10)
time = linspace(-1.4e-6, 6.4e-6, 10)
plot([0.260e-6]*len(amp), amp, '--', lw = 1.5, color = 'grey')
plot([0.520e-6]*len(amp), amp, '--', lw = 1.5, color = 'grey')
plot([2e-6]*len(amp), amp, '--', lw = 1.5, color = 'grey')
plot(time , [1.43]*len(time), '--', lw = 1.5, color = 'grey')
"""

grid(True, which = 'both')
plot(tDDS, rDDS, color = 'grey', lw = 1.5, label = "$C_{reference}(t)\ (Trigger)$")
plot(taom0, raom0, lw = 1.5, label = r"$C_{Optics}(t)$")
plot(tVCO, rVCO, lw = 1.5, label = r"$C_{Electronics}(t)$")

#plot(tVCO, rVCO)
#plot(taom1, raom1)
ylim(0, 1.5)
xlim(-1.4e-6, 6.4e-6)
p.scale_xaxis(ax,1e6)
xlabel(r"Time $\mu$s")
ylabel(r"Amplitude relative to Set Value")
yticks(arange(0.0,1.6,0.2), [r'0%', '20%', '40%', '60%', '80%', r'100%', '120%', '140%'])

legend(loc = 4)
#grid()
tight_layout()
savefig("Final_Step_Response_Electronics.png")
show() 

# <codecell>

### Derivative of Step Response ###

figure()
ax = subplot(111)
p.grey_title("Derivative of the Step Response Function")
p.colorcycle(3)

Taom0, Daom0 = dy_dx(taom0, raom0)
plot(Taom0, Daom0, lw = 1.5, label = r'$\frac{d}{dt}C_{Optics}(t)$')

Taom1, Daom1 = dy_dx(taom1, raom1)
#plot(Taom1, Daom1, lw = 1.5, label = 'AOM1 Derivative')

TVCO, DVCO = dy_dx(tVCO, rVCO)
plot(TVCO, DVCO, lw = 1.5, label = r'$\frac{d}{dt}C_{Electronics}(t)$')

xlim(0, 6.4e-6)
p.scale_xaxis(ax,1e6)
xlabel(r"Time $\mu$s")
ylabel(r"Amplitude relative to Set Value")

legend(loc = 4)
grid()
tight_layout()
savefig("Final_Step_Response_Derivative.png")
show()

# <markdowncell>

# The system transfer function $H(s)$ is obtained by Fourier transformation of the derivative
# $$ H(2\pi i f) = \mathcal F \lbrace \frac{d}{dt} C(t) \rbrace.$$
# 
# The -3 dB value from the system transfer function indicates the bandwidth of the loop. It can be found at 
# $$ f_{electronics} ^{-3dB}= 1.35\ \textrm{MHz} \qquad f^{-3dB}_{optic} = 1.2\ \textrm{MHz}$$
# 
# and is therefore not very different to the one obtained for the electronics alone. The largest influence of the AOM consists in the phase decrease by the dead time. At 1 MHz the phase of the system transfer function with the AOM is decreased by 70$^\circ$ compared to the one where only the electronics with negligible dead times contribute.

# <codecell>

### FFT(Derivative) ###
    
faom0, FFTaom0 = frequency_analysis(Taom0, Daom0)
faom1, FFTaom1 = frequency_analysis(Taom1, Daom1)
fVCO, FFTVCO = frequency_analysis(TVCO, DVCO)

SystemTransfer0 = (FFTaom0*1/abs(FFTaom0[0]))
SystemTransfer1 = (FFTaom1*1/abs(FFTaom1[0]))
SystemTransferVCO = (FFTVCO*1/abs(FFTVCO[0]))

xmin, xmax = 8.5e3, 2.6e6#)#90e3,2.8e6)
ymin, ymax = -15,5

figure()
p.grey_title(r"System Transfer Function $H(s)$")
p.colorcycle(3)
Plot_abs = subplot(111)

freqs = linspace(xmin,xmax,100)
amps = linspace(ymin, ymax, 100)
plot(freqs, [-3]*len(freqs), '--', lw = 2, color = 'grey')
plot(freqs, [0]*len(freqs), '--', lw = 2, color = 'grey')
plot([1.2e6]*len(amps), amps, '--', lw = 2, color = 'grey')
plot([1.35e6]*len(amps), amps, '--', lw = 2, color = 'grey')

plot(list(faom0), list( 10*log10(abs(SystemTransfer0))), lw = 1.5, label = r'$|H_{Optics}|$')
#plot(faom1, 10*log10(abs(SystemTransfer1)), lw = 1.5, label = "FFT(Derivative AOM1) ")
plot(fVCO, 10*log10(abs(SystemTransferVCO)), lw = 1.5, label = r'$|H_{Electronics}|$')
    
xscale("log")
xlabel("Frequency $f$ [Hz]")
grid(True, which="both")
ylabel("Amplitude [dB]")

legend(loc=3)
ylim(ymin, ymax)

Plot_phase = Plot_abs.twinx()

p.colorcycle(3)

import numpy as np
Phase0 = np.unwrap(np.angle(SystemTransfer0))*180/pi #np.angle(FFT[index], deg = True)#
plot(faom0, array(Phase0),'--', lw = 1.5, label = r'$arg(H_{Optics})$')

Phase1 = np.unwrap(np.angle(SystemTransfer1))*180/pi #np.angle(FFT[index], deg = True)#
#plot(faom1, Phase1,'--', lw = 1.5, label = "Phase 1")

PhaseVCO = np.unwrap(np.angle(SystemTransferVCO))*180/pi #np.angle(FFT[index], deg = True)#
plot(fVCO, array(PhaseVCO),'--', lw = 1.5, label = r'$arg(H_{Electronics})$')
    
ylim(-220, 20)
ylabel(r"Phase in $^\circ$")
legend(loc=2)
xlim(xmin, xmax)
#tight_layout()
savefig("Final_System_Transfer_Functions_Opt_and_Electronics.png")
show()

# <codecell>

OpenLoop = lambda ClosedLoop: ClosedLoop/(1-ClosedLoop)
    
OpenLoop0 = OpenLoop(SystemTransfer0)
OpenLoop1 = OpenLoop(SystemTransfer1)
OpenLoopVCO = OpenLoop(SystemTransferVCO)

# <codecell>

wn = 7.4e6
D = 0.27
H = lambda s, D: wn**2/(s**2+2*D*wn*s+wn**2)
    
OpenLoopAn = OpenLoop(H(2*pi*1j*fVCO, D))

# <markdowncell>

# From the system transfer function $H(s)$ the open loop gain can be calculated
# $$ G(s) = \frac{H(s)}{1-H(s)}.$$
# 
# Its absolute value and phase are plotted in figure \ref{fig:FIG}. 
# 
# The crossing frequency for which the open loop gain $|G(s)|$ is unity, lies at 487 kHz (INTERPRETATION FOR THIS?). The stability of the system depends on the phase margin at this point, which can be read off from the plot as
# 
# $$ \varphi_{electronics}^{margin} = 180^\circ -110^\circ = (70\pm 2)^\circ \qquad \varphi^{margin}_{optics} = 180^\circ - 145^\circ = (35 \pm 2)^\circ .$$
# 
# The value for the phase margin of the optics estimated by the overshoot in the step response was $(27\pm 2)^\circ$, which is not in agreement, but very close to the actual value. 
# 
# The additional dead time caused by the AOM was measured in section \ref{REF} to be $T_{AOM,\ dead} = 235\ ns$. At the cross frequency this should cause an additional phase compared to the electronic circuit of
# $$ \Delta \varphi_{expected}^{additional}= 2\pi \cdot 487\textrm{ kHz} \cdot 235\textrm{ ns} \approx 40^\circ \qquad  \Delta \varphi^{additional}_{measured} = 145^\circ - 110^\circ = 35^\circ ,$$
# 
# which shows that the additional phase is in the expected magnitude. 

# <codecell>

### Open Loop Transfer Function ###

figure()
p.grey_title(r"Open Loop Transfer Function $G(s)$")
p.colorcycle(3)
Plot_abs = subplot(111)
ymin, ymax = -15,15

amps = linspace(ymin, ymax,10)
plot(freqs, [0]*len(freqs), '--', lw = 2, color = 'grey')
plot([487e3]*len(amps), amps, '--', lw = 2, color = 'grey')

plot(faom0, 10*log10(abs(OpenLoop0)), lw = 1.5, label = r'$|G_{Optics}|$')
#plot(faom1, 10*log10(abs(OpenLoop1)), lw = 1.5, label = "FFT(Derivative AOM1) ")
plot(fVCO, 10*log10(abs(OpenLoopVCO)), lw = 1.5, label = r'$|G_{Electronics}|$')
#plot(fVCO, 10*log10(abs(OpenLoopAn)), lw = 1.5, color = 'black', label = "FFT(Derivative An) ")

xscale("log")
xlabel("Frequency in Hz")
grid(True, which="both")
ylabel("Amplitude in dB")

legend(loc=3)
ylim(ymin,ymax)

Plot_phase = Plot_abs.twinx()

p.colorcycle(3)

Phase0 = np.unwrap(np.angle(OpenLoop0))*180/pi
#plot(faom0, Phase0,'--', lw = 1.5, label = "Phase 0")
#Phase1 = np.unwrap(np.angle(OpenLoop1))*180/pi
f,O = data_region(faom1, OpenLoop1,140e3,5e6)#.5e5,5e6)
P = np.unwrap(np.angle(O))*180/pi
plot(f, P,'--', lw = 1.5, label = r'$arg(G_{Optics})$')

plot(freqs, [-180]*len(freqs), '--', lw = 2, color = 'grey')
plot(freqs, [-145]*len(freqs), '--', lw = 1.5, color = 'grey')
plot(freqs, [-110]*len(freqs), '--', lw = 1.5, color = 'grey')

PhaseVCO = np.unwrap(np.angle(OpenLoopVCO))*180/pi
plot(fVCO, PhaseVCO, '--', lw = 1.5, label = r'$arg(G_{Electronics})$')

PhaseAn = np.unwrap(np.angle(OpenLoopAn))*180/pi
#plot(fVCO, PhaseAn, '--', lw = 1.5, color = 'black', label = "Phase An")
    
ylim(-220, -100)
ylabel(r"Phase in $^\circ$")
legend(loc=4)
xlim(xmin,xmax)
#tight_layout()
savefig("Final_Open_Loop_Transfer_Function_Opt_and_Electronics.png")
show()

# <codecell>


