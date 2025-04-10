import matplotlib.pyplot as plt
import sys
#import numpy as np
plt.style.use('seaborn-whitegrid')

ax = plt.axes()

xname=[]
yname=[]

x = []
y = []

opt =1

if opt == 0:
  pl_col_0 = 1
  pl_col_1 = 2
elif opt == 1:
  pl_col_0 = 8
  pl_col_1 = 9
else:
  pl_col_0 = 5
  pl_col_1 = 6 

col = 0

#for line in open ('../HPCresults/IbarRes3Nstag50ch_3_dis6/stress_strain.txt', 'r'):
#for line in open ('../results/stress_strain.txt', 'r'):
for line in open (sys.argv[1], 'r'):
#for line in open ('../stress_strain.txt', 'r'):
#for line in open ('../HPCresults/Ibardamagedis7/results/stress_strain.txt', 'r'):
  lines = [i for i in line.split("\t")]
  if (col == 0):
    xname.append(lines[pl_col_0])
    yname.append(lines[pl_col_1])
    col = 1
    continue
  x.append(float(lines[pl_col_0]))
  y.append(float(lines[pl_col_1]))
  
plt.title("Stress Strain curve")
plt.xlabel(xname[0])
plt.ylabel(yname[0])

#plt.plot(x[1:-1], y[1:-1], ':r')
plt.plot(x[0:-1], y[0:-1],'o:r',ms=4)
plt.show()
