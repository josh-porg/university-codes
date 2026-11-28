"""DMD book chapter 6 (``Algorithm_6_1``, ``Algorithm_Sec_6_2``): DMD with control.

An unstable system ``x_{k+1} = A x_k + B u_k`` with A = diag(1.5, 0.1),
B = (1, 0) and proportional feedback u = -x1 is identified from 20 steps.
Plain DMD returns the closed-loop matrix (stable, wrong); DMDc with known
B recovers A; DMDc with unknown B (section 6.2, :func:`unicodes.decomposition.dmdc`)
identifies both from the state and input data. (The book's section 6.2
code took only one row of the input block of the singular vectors, which
is wrong for more than one input; the library uses all of them.)
"""

import numpy as np

from unicodes.decomposition import dmdc

A = np.array([[1.5, 0], [0, 0.1]])
B = np.array([[1.0], [0]])
K = -1.0
m = 20
X = np.zeros((2, m + 1))
X[:, 0] = [4, 7]
U = np.zeros((1, m))
for j in range(m):
    U[:, j] = K * X[0, j]
    X[:, j + 1] = A @ X[:, j] + (B * U[:, j]).ravel()

Xk, Xp = X[:, :-1], X[:, 1:]
Uu, S, Vh = np.linalg.svd(Xk, full_matrices=False)
r = int(np.sum(S > 1e-10))
Uu, S, V = Uu[:, :r], S[:r], Vh[:r].T
print("A_DMD (sees the closed loop A + B K):\n", Xp @ V / S @ Uu.T)
print("A_DMDc with B known:\n", (Xp - B @ U) @ V / S @ Uu.T)

# B unknown: with u = K x1 the inputs are collinear with the states, so the
# data cannot separate A from B; add a small excitation as in the book's section 6.2 setting.
rng = np.random.default_rng(0)
X2 = np.zeros((2, m + 1))
X2[:, 0] = [4, 7]
U2 = np.zeros((1, m + 1))
for j in range(m):
    U2[:, j] = K * X2[0, j] + rng.standard_normal()
    X2[:, j + 1] = A @ X2[:, j] + (B * U2[:, j]).ravel()
res = dmdc(X2, U2, tol=1e-10)
print("DMDc with B unknown (excited input): A =\n", res.A, "\nB =\n", res.B)
