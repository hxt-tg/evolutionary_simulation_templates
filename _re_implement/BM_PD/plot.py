import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple

Conf = namedtuple("Conf", "file color")

configures = [
    Conf("output/p/beta=0.20_A=0.50.csv", 'b'),
    Conf("output/p/beta=0.20_A=1.50.csv", 'g'),
]

epsilon = 0.2

for c in configures:
    p = np.loadtxt(c.file, skiprows=1, delimiter=",", usecols=(2, 3))
    p = np.average(p.reshape((-1, 100, 2))[:, :, 0], axis=0)
    p = (1-p)*epsilon + (1-epsilon)*p

    plt.plot(p, c.color)

plt.xlim([0, 100])
plt.ylim([0, 0.6])
plt.show()