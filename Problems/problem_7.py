import numpy as np
import matplotlib.pyplot as plt
import pyarma


filename = "solution.bin"
solution = pyarma.mat()
solution = np.array(solution.load(filename))