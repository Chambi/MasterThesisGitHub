# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Joined Function and RBW Value Readout

# <headingcell level=2>

# Load Traces into Arrays, Plot and Read out the Resolution Bandwidth RBW

# <codecell>

my_name			= 'Joined_Function.csv'

import numpy as np # numpy for array handling and numeric calculations
import panna as p # panna will search for all csv files and read them into one data array
from pylab import *

data = p.load_all()
len_files = p.number_of_files()
RBW = p.load_RBW()

# <codecell>

x, y = [0]*len_files, [0]*len_files
for i in arange(0,len_files):
    x[i], y[i] = array(data[i][:,0]), array(data[i][:,1])

# <codecell>

def center(x_values, y_values):
	y_values = list(y_values)
	maximum_index 	= y_values.index(max(y_values))	
	centerfreq	= x_values[maximum_index]
	return centerfreq

# <codecell>

figure()
subplot(121)
p.grey_title("Data on Linear Scale")
p.colorcycle(len_files)
for i in arange(0, len_files):
    plot(x[i], y[i], label = str(i))
p.scale_xaxis(subplot(121),1e-6)
xlabel(r'Frequency [MHz]')
ylabel(r'Amplitude [dB]')

subplot(122)
p.grey_title("Data on Log Scale")
p.colorcycle(len_files)
for i in arange(0, len_files):
    plot(x[i]-center(x[i],y[i]), y[i], label = 'RBW='+str(RBW[i]))
xscale("log")
xlabel(r'Offset Frequency [Hz]')
ylabel(r'Amplitude [dB]')
leg = legend()
leg.get_frame().set_alpha(0.35)

tight_layout()
savefig("Step_1_Plot_all_Data_in_Folder.png")
show()

# <headingcell level=2>

# Normalize to 1 Hz Resolution Bandwidth

# <codecell>

for i in arange(0, len_files):
    y[i] = array(y[i])-10*log10(RBW[i]) 

# <codecell>

figure()

subplot(121)
p.grey_title("Data on Linear Scale")
p.colorcycle(len_files)
for i in arange(0, len_files):
    plot(x[i], y[i], label = str(i))
ylabel("Amplitude [dBm/Hz]")
xlabel("Frequency [Hz]")

subplot(122)
p.grey_title("Data on Log Scale")
p.colorcycle(len_files)
for i in arange(0, len_files):
    plot(x[i]-center(x[i],y[i]), y[i], label = 'RBW='+str(RBW[i]))
xscale("log")
xlabel(r'Offset Frequency [Hz]')
ylabel(r'Amplitude [dBm/Hz]')
leg = legend()
leg.get_frame().set_alpha(0.35)
tight_layout()
savefig("Step_2_RBW_normalized_Traces.png")
show()

# <headingcell level=2>

# Joined Function

# <codecell>

############################ Finding the Joined Function ################################
def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx

# <codecell>

minimal, maximal = [], []
for i in arange(0, len_files):
    minimal.append([i, x[i][0]])
    maximal.append([i, x[i][-1]])

minimal = array(sorted(minimal, key=lambda element: element[1]))
maximal = array(sorted(maximal, key = lambda element: element[1]))
print minimal[:,0].astype(int), maximal[:,0].astype(int)

# <codecell>

joined_x = []
joined_y = []

for i in minimal[:,0].astype(int):
    if i == minimal[:,0][-1]:
        jumping_index = len(x[i])
    else: 
        jumping_index = find_nearest(x[i],minimal[i+1,1]) 
    for j in arange(0, jumping_index):
        joined_x.append(x[i][j])
        joined_y.append(y[i][j])
    print "appending array ",i," from 0 till ", jumping_index

for i in maximal[:,0].astype(int):
    jumping_index = find_nearest(x[i],joined_x[-1])
    for j in arange(jumping_index, len(x[i])):
        joined_x.append(x[i][j])
        joined_y.append(y[i][j])
    print "appending array ",i," from ", jumping_index, " till ", len(x[i])
##########################################################################################################

# <codecell>

figure()

subplot(121)
p.grey_title("Data on Linear Scale")
p.colorcycle(len_files)
for i in arange(0, len_files):
    plot(x[i], y[i], label = str(i))
    
plot(joined_x, joined_y, lw = 0.5, color = 'black')
ylabel("Amplitude [dBm/Hz]")
xlabel("Frequency [Hz]")

subplot(122)
p.grey_title("Data on Log Scale")
p.colorcycle(len_files)
for i in arange(0, len_files):
    plot(x[i]-center(x[i],y[i]), y[i], label = 'RBW='+str(RBW[i]))
    
xlabel(r'Offset Frequency [Hz]')
ylabel(r'Amplitude [dBm/Hz]')
leg = legend()
leg.get_frame().set_alpha(0.35)

plot(joined_x-center(x[0],y[0]), joined_y, lw = 0.5, color = 'black')
xscale("log")
tight_layout()
savefig("Step_3_Finding_Joined_Function.png")
show()

# <headingcell level=2>

# Output Joined Function into File

# <codecell>

filename_my_joined = my_name
FILE_my_joined = open(filename_my_joined, "w")

my_try = []
for i in range(0,len(joined_x)):
	my_try.append(str(joined_x[i])+","+str(joined_y[i]))

for item in my_try:
	print >> FILE_my_joined, item
    
print "Joined Function saved in ", my_name

# <codecell>


