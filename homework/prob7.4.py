from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 2:
    print('usage: python prob7.4.py n')

n = int(sys.argv[1])

def lane_emden(x, Y):
    y = Y[1]
    z = -Y[0]**n - (2/x)*Y[1]
    return y,z

end_xi = 15
points = np.linspace(1, 15, 10000)
r = solve_ivp(lane_emden, (1, end_xi), (1, 0), t_eval = points)

plt.style.use('seaborn')
linestyles = ['-', '--']
plt.figure(figsize=(10,16))
[plt.plot(r.t, y, linestyle=ls) for (y, ls) in zip(r.y, linestyles)]
labels = [r'y', r'z']
plt.legend(labels, loc=1)
plt.xlabel('x')
plt.ylabel('y and z')
plt.show()