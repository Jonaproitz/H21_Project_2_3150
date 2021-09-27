import numpy as np
import matplotlib.pyplot as plt
import pyarma


filename_1 = "solution.bin"
solution = pyarma.mat()
solution.load(filename_1)

sol = np.transpose(np.array(solution))
x = sol[0]
V = sol[1:]

filename_2 = "analytic_solution.bin"
analytic_solution = pyarma.mat()
analytic_solution.load(filename_2)

U = np.transpose(np.array(analytic_solution))

plt.figure()

for i in range(3):
    v = V[i]
    u = U[i]
    if np.allclose(v, -u):
        v = -v
    plt.plot(x, v)
    plt.plot(x, u, "--")

plt.savefig("Eigenvec.pdf")