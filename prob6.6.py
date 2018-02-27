from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np


def abundances_plus_30(t, y):
  return abundances(t, y, 1.3)

def abundances_minus_30(t, y):
  return abundances(t, y, .7)

def abundances(t, y, alpha_12 = 1):
  rho = 1e4
  Na = 6.022e23
  A4 = 4
  A12 = 12
  A16 = 16
  A20 = 20
  a3 = 1.48e-18 / Na**2
  a12 = 3.03e-17 * alpha_12 / Na
  a16 = 1.05e-22 / Na
  sec_per_yr = 365.25 * 24 * 3600

  x4, x12, x16, x20 = y

  RHS_sec = [
    rho * Na * x4 * (-.5 * rho * Na * x4**2 * a3 / A4**2 - x12*a12/A12 - x16 * a16 / A16),
    rho * Na / A4 * x4 * ((1/6) * rho * Na * A12 / A4**2 * x4**2 * a3 - x12 * a12),
    rho * Na / A4 * x4 * (x12 * a12 * A16 / A12 - x16 * a16),
    rho * Na / A4 * x4 * A20 * x16 * a16 / A16
  ]

  RHS_yr = [x * sec_per_yr for x in RHS_sec]
  return RHS_yr

end_time = 1e7
times = np.linspace(0, end_time, 5000)
res = solve_ivp(abundances, (0, end_time), (1,0,0,0), t_eval = times)
res_plus_30 = solve_ivp(abundances_plus_30, (0, end_time), (1,0,0,0), t_eval = times)
res_minus_30 = solve_ivp(abundances_minus_30, (0, end_time), (1,0,0,0), t_eval = times)

plt.style.use('seaborn')
linestyles = ['-', '--', '-.', ':']
plt.subplot(311)
[plt.plot(res.t, y, linestyle=ls) for (y, ls) in zip(res.y, linestyles)]
labels = [r'$X_4$', r'$X_{12}$', r'$X_{16}$', r'$X_{20}$']
plt.legend(labels, loc=1)
plt.title(r'$\langle\alpha 12 \rangle$ = 100%')

plt.subplot(312)
[plt.plot(res_plus_30.t, y, linestyle=ls) for (y, ls) in zip(res.y, linestyles)]
labels = [r'$X_4$', r'$X_{12}$', r'$X_{16}$', r'$X_{20}$']
plt.legend(labels, loc=1)
plt.title(r'$\langle\alpha 12 \rangle$ = 130%')

plt.subplot(313)
[plt.plot(res_minus_30.t, y, linestyle=ls) for (y, ls) in zip(res.y, linestyles)]
labels = [r'$X_4$', r'$X_{12}$', r'$X_{16}$', r'$X_{20}$']
plt.legend(labels, loc=1)
plt.title(r'$\langle\alpha 12 \rangle$ = 70%')





plt.suptitle(r'Abundances from $3\alpha$', fontsize=14)
plt.subplots_adjust(left=0.05, bottom=0.05,right=0.99,top=0.95,wspace=0.14,hspace=0.14)
plt.show()
