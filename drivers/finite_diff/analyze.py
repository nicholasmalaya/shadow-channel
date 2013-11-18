import os
from pylab import *
from numpy import *
from numpy import polynomial

u_plus, u_minus = [], []

for i in range(1, 1000000):
    f_plus = 'Re2100_{0:06d}.txt'.format(i)
    f_minus = 'Re1900_{0:06d}.txt'.format(i)
    if not (os.path.exists(f_plus) and os.path.exists(f_minus)):
        break
    print('reading ', f_plus, f_minus)
    up, um = loadtxt(f_plus), loadtxt(f_minus)
    u_plus.append(0.5 * (up + up[::-1]))
    u_minus.append(0.5 * (um + um[::-1]))

um_plus, um_minus = array(u_plus).mean(0), array(u_minus).mean(0)
ustd_plus, ustd_minus = array(u_plus).std(0) / sqrt(len(u_plus)), \
                        array(u_minus).std(0) / sqrt(len(u_minus))
du = (um_plus - um_minus) / 200.
du_std = sqrt(ustd_plus**2 + ustd_minus**2) / 200.

y, w = polynomial.legendre.leggauss(du.shape[0])

figure(figsize=(12, 8))
titles = ['U', 'V', 'W', 'U2', 'V2', 'W2']
for i in range(du.shape[1]):
    subplot(2,3,i+1)
    plot(y, du[:,i], '-k')
    plot(y, du[:,i] + du_std[:,i] * 2, '--k')
    plot(y, du[:,i] - du_std[:,i] * 2, '--k')
    title(titles[i])

savefig('fintie_diff')
