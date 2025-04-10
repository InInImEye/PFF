# To solve Ax = b
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

# Reads RHS
def read_rhs(file_path):
  rhs = np.loadtxt(file_path, dtype=float)
  
  return rhs

# Reads Matrix
def read_sparse_matrix(file_path):
    rows = []
    cols = []
    values = []

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            # Parse the row, col and value
            row, col = map(int, parts[0][1:-1].split(','))
            value = float(parts[1])
            # Store the row, col and value
            rows.append(row)
            cols.append(col)
            values.append(value)

    # Construct the COO sparse matrix
    sparse_matrix = coo_matrix((values, (rows, cols)))

    return sparse_matrix

b = read_rhs("System_rhs.txt")
A = read_sparse_matrix("System_matrix.txt")
x = spsolve(A, b)
np.allclose(A.dot(x).toarray(), b.toarray())

print(x)
