import numpy as np
import matplotlib.pyplot as plt 
import os
import scipy.optimize


N = [N for N in range(2, 61)]
iterations = []

for N_i in N:
    os.system(f"./Problem_6.exe {N_i}")
    with open("iterations.txt") as file:
        iterations.append(int(file.readline()))

N = np.array(N)
iterations = np.array(iterations)

def curve(N, a, b, c):
    return a*N**2 + b*N + c 

popt, pcov = scipy.optimize.curve_fit(curve, N, iterations)

print(popt)

a, b, c = popt


plt.figure()
plt.title("Change in number of iterations as a function of N", fontsize=15)

plt.plot(N, iterations, label="Number of iterations")
plt.plot(N, curve(N, a, b, c), label="Fitted curve")

plt.xlabel("N", fontsize=12)
plt.ylabel("Iterations", fontsize=12)
plt.legend()
plt.savefig("../LaTeX/Iterations.pdf")