#Python file for plotting results
import os 

from pylab import *
with open('results.txt', 'r') as f:
	line = f.readline()
	data0 = line.split()
	number_of_lines = int(data0[0])
	probability = zeros([3,number_of_lines])
	rho = zeros(number_of_lines)
	line = f.readline()
	eigen_value = double(line)
	f.readline()
	f.readline()
	for i in range(0,number_of_lines):
		line = f.readline()
		data0 = line.split()
		rho[i] = data0[3]
		probability[0][i] = data0[0]
		probability[1][i] = data0[1]
		probability[2][i] = data0[2]


plot(rho,probability[0][:], rho,probability[1][:], rho,probability[2][:])

xlabel(r'$\rho$', fontsize = 20)
ylabel('Probability distribution', fontsize = 20)
xticks(fontsize = 14) 
yticks(fontsize = 14)
legend(['Gound state', 'First exitation', 'Second exitation'], fontsize = 16)
filename = "omega5.eps"
savefig(filename)
os_string = "mv " + filename + " ~/Skole/H15/FYS3150/Projects/Project2/latex/."
print os_string
os.system(os_string)
show()



