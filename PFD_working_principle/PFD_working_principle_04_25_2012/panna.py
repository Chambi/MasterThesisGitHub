"""
import panna as p
from pylab import *

#p.plot_all()
data = p.load_all()
p.show_n()
len_files = p.number_of_files()
x,y,z,a,b = [0]*len_files, [0]*len_files, [0]*len_files, [0]*len_files, [0]*len_files

for i in arange(0, len_files):
    len_columns = p.number_of_columns(i)
    for j in arange(0, len_columns-1):
        if j == 0: x[i] = array(data[i][:,j])
        if j == 1: y[i] = array(data[i][:,j])
        if j == 2: z[i] = array(data[i][:,j])
        if j == 3: a[i] = array(data[i][:,j])
        if j == 4: b[i] = array(data[i][:,j])
"""

import glob
from pylab import *

# access all csv files:
files = glob.glob('*.csv')
# a line which contains data, only contains this signs:
data_signs = ['\n','',' ','.',',','0','1','2','3','4','5','6','7','8','9','e','E','-','+']

# lines[i] oontains the whole content of file files[i]
def load_lines(file):
	lines = []
	f = open(file,'r')
	lines.append(f.readlines()[0:])
	#print "file ",file," found "
	#print "loaded lines ", lines[0][0:3], "..."
	f.close()
	return lines[0]

# The Spectrum Analyzer saves the Resolution Bandwidth into the header. It is to be read out here:
def get_RBW(datafile):
	read = load_lines(datafile)
	line_array = [line.split(',') for line in read]
	RBW_element = float(line_array[11][1]) # SA saves RBW there
	return RBW_element

def number_of_files():
	return len(files)

def load_data(i):
	#print files[i]," access data via index ",i 
	lines = load_lines(files[i])
	data = []
	first_line = True
	for line in lines[0:-2]:
	# if line contains anything but numbers . and , skip, otherwise append
		if all([element in data_signs for element in line]):
			if first_line:
				#print "The first line of file ",i," is ",str(line)
				first_line = False
			line_array = line.strip('\n')
			line_array = line.split(',')
			if ('' not in line_array):
				for element in line_array:
					line_array = [float(element) for element in line_array]
					data.append(line_array)	
	return array(data)
	
def load_all():
	i = 0
	data = [0]*len(files)
	for file in files:
		data[i] = load_data(i)
		i += 1
	return data

def load_RBW():
	i = 0
	RBW = [0]*len(files)
	for file in files:
		RBW[i] = get_RBW(file)
		i+=1
	return RBW
	
def colorcycle(num_plots):
	colormap = plt.cm.gist_rainbow
	plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.35, 1.0, num_plots)])

def grey_title(my_title):
	title(my_title, style='italic', color = 'white', bbox={'facecolor':'grey', 'alpha':0.5})

def plot_all():
	data = load_all()
	i = 0
	for file in files:
		figure()
		columns = shape(data[i])[1]
		print "file ",file," has ",columns," columns"
		colorcycle(columns)
		grey_title(str(files[i]))
		for j in arange(1,columns):
			print "plotting column ",j
			plot(data[i][:,0], data[i][:,j], label = "column "+str(j))
		legend()
		savefig("plot_data_"+str(files[i])+".png")
		close()
		i += 1
		
def show_n():
	i = 0
	for file in files:
		print "access file ",file," via index ",i
		i += 1
	
def number_of_columns(i):
	data = load_data(i)
	return shape(data)[1]
	
def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx

def normalize(array):
	minimum = min(array)
	new_array = array - minimum
	maximum = max(new_array)
	result = new_array/maximum
	return result
	
	

