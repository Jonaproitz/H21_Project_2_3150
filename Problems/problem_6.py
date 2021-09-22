import numpy as np
import matplotlib.pyplot as plt 
import os


N = [N for N in range(2, 10)]
iterations = []

for N_i in N:
    os.system(f"~/Documents/Sem_7/FYS3150/Project_2/Problems/Problem_6.exe {N_i}")
    with open("iterations.txt") as file:
        iterations.append(int(file.readline()))

