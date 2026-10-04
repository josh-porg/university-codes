"""Math 290 and Math 590 linear algebra (``Linear_Algebra_Experimentation``, ``find_basis``,
``linear_algebra_test_2``, ``linear_algebra_test_final``, ``linear_algebra_hw_7``), with sympy.
"""

from itertools import combinations

import numpy as np
import sympy as sp

# Experimentation: ||A|| = ||A^T|| and the SVD of a 3x2 matrix from eig(A^T A), eig(A A^T)
A = np.random.default_rng(0).random((5, 99))
print("||A||_2 - ||A^T||_2 =", np.linalg.norm(A, 2) - np.linalg.norm(A.T, 2))
A = sp.Matrix([[1, -2], [-3, 0], [1, 2]])
print("A^T A =", A.T * A, " eigenvalues of A A^T:", (A * A.T).eigenvals(), " singular values:", A.singular_values())

# find_basis: which 4 of the 6 vectors form a basis of R^4 (non-zero determinant)
V = sp.Matrix([[3, 3, 1, -3], [3, 6, 3, -2], [3, 6, 3, -1], [-3, -6, -1, 3], [-3, -6, -3, 5], [0, -3, 0, -1]]).T
print("rank of the 6 vectors:", V.rank())
for idx in combinations(range(6), 4):
    if V[:, list(idx)].det() == 0:
        print(f"  vectors {[i + 1 for i in idx]} are dependent")

# Test 2: eigenvectors and the span of two of them
A = sp.Matrix([[-3, 18, 6], [0, -9, -2], [0, 24, 5]])
P, D = A.diagonalize()
print("Test 2: eigenvalues", list(D.diagonal()), " eigenvectors", [list(P[:, i]) for i in range(3)])
print("        column space of eigenvectors 2 and 3:", [list(v) for v in P[:, 1:].columnspace()])

# Final: least squares
A = sp.Matrix([[3, -1], [-1, 2], [0, -3]])
b = sp.Matrix([1, -3, 2])
print("Final: least-squares solution", list((A.T * A).inv() * A.T * b))

# Math 590 HW 7
A = sp.Matrix([[0, -1], [1, 0]])
print("590 HW 7: eig", A.eigenvals(), " LU", A.LUdecomposition()[:2], " QR", A.QRdecomposition())
x = sp.Symbol("x")
print("  roots of x^2 - 7x + 10 + 18i:", sp.solve(x**2 - 7 * x + 10 + 18 * sp.I, x))
B = sp.Matrix([[2, -3 - 3 * sp.I], [3 + 3 * sp.I, 5]])
print("  B normal:", sp.simplify(B * B.H - B.H * B) == sp.zeros(2), " Hermitian:", B == B.H, " eig:", B.eigenvals())
x1, x2, x3, x4 = sp.symbols("x1:5")
T = sp.Matrix([[x1, x2], [x3, x4]])
eqs = list(T - T.T) + list(T**2 - T) + list(T**2 - T.T) + list((T.T) ** 2 - T)
print("  symmetric idempotent T:", sp.solve([e for e in eqs if e != 0], [x1, x2, x3, x4], dict=True))
