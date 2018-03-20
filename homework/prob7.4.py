from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import sys

n = 2

def lane_emden(x, Y):
    global n
    y, z = Y

    dydx = z
    if x == 0:
        dzdx = -y**n + (2/3)
    else:
        dzdx = -y**n -  (2/x)*z

    return dydx, dzdx

end_xi = 15
points = np.linspace(0, end_xi, 1000)



plt.style.use('seaborn')
linestyles = ['-', '--']
plt.figure(figsize=(10,16))
for n in np.linspace(2, 4, 7):
    r = solve_ivp(lane_emden, (0, end_xi), (1, 0), t_eval = points)
    mask = (r.y[0] >= 0) & (r.y[1] <= 1)
    x = r.t[mask]
    y = r.y[0][mask]
    z = r.y[1][mask]
    plt.subplot(211)
    plt.plot(x, y, label=f'n={n:.1f}')
    plt.subplot(212)
    plt.plot(x, z, label=f'n={n:.1f}')
plt.subplot(211)
plt.xlabel(r'$x=\xi$')
plt.ylabel(r'$y=\theta_n$')
plt.subplot(212)
plt.xlabel(r'$x=\xi$')
plt.ylabel(r'$z=\frac{d\theta_n}{d\xi}$')
plt.legend(loc = 4)
plt.suptitle(r'y and z vs x for a range of n $\in$ [2,4]')
plt.savefig('prob7.4.pdf', format='pdf')