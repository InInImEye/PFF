# Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html#scipy.sparse.linalg.spsolve

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
B = csc_matrix([[2, 0],  [-1, 0], [2, 0]], dtype=float)
x = spsolve(A, B)
np.allclose(A.dot(x).toarray(), B.toarray())
print(x)
