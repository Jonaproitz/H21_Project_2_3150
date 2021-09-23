import numpy as np
import matplotlib.pyplot as plt
import pyarma


filename = "solution.bin"
solution = pyarma.mat()
solution.load(filename)

sol = np.transpose(np.array(solution))
x = sol[0]
V = sol[1:]

