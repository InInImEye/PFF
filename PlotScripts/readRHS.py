import numpy as np

def read_rhs(file_path):
  rhs = np.loadtxt(file_path, dtype=float)
  
  return rhs

b = read_rhs("System_rhs.txt")
print(b)
print(np.size(b))
