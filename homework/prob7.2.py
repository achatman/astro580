output = '  M (M_solar),  R (R_solar),Rho_c (g/cm3),<Rho> (g/cm3), Rho_c/<Rho> \n'
with open('prob7.2.csv', mode='r') as f:
  for line in f:
    M,R,rhoc = line.strip().split(',')
    Mx = float(M) * 1.989e33
    Rx = float(R) * 1e10
    rhoc = float(rhoc)
    ave_rho = 3 * Mx / (4 * 3.14 * Rx**3)
    outstr = f'{M:>13},{R:>13},{rhoc:13.3},{ave_rho:13.3},{rhoc/ave_rho:13.3}'
    output += outstr + '\n'

with open('prob7.2.txt', mode='w') as f:
  f.write(output)
