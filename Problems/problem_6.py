import numpy as np
import matplotlib.pyplot as plt 
import os
import scipy.optimize


N = [N for N in range(2, 61)]
iterations = []

for N_i in N:
    os.system(f"~/Documents/Sem_7/FYS3150/Project_2/Problems/Problem_6.exe {N_i}")
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

plt.plot(N, iterations)
plt.plot(N, curve(N, a, b, c))

plt.show()