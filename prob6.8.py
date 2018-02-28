from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

def temp(t, y):
  T = y[0]
  T9 = T * 10e-9
  rho = 2e5
  Y = 1
  mu_e = 2#9.11e-27
  #me = 2.74e-3
  B = 9.739e5
  x = (rho / (mu_e * B))**.33
  cve = 1.35e5 / rho * T * x * (1 + x**2)**.5
  cvHe = 9.35e7#4.125e7
  cv = cvHe + cve
  ep3a = 5.1e8 * rho**2 * Y**3 / T9**3 * np.exp(-4.4027 / T9)

  dTdt = ep3a / (cve + cv)

  return dTdt * 86400

end_time = 1e11
low_times = np.linspace(0, 20, 10000, endpoint=False)
times = low_times#np.append(low_times, np.linspace(20, end_time, 10000))
res = solve_ivp(temp, (0, end_time), (1.5e8,), t_eval = times)
plt.style.use('seaborn')
plt.figure(figsize=(10,16))
plt.plot(res.t, res.y[0])

plt.show()



