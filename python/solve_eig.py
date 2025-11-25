import numpy as np
import scipy.io as sio
import scipy.linalg as sla

import matplotlib.pyplot as plt

import mpmath as mp
from python.gen_eig_ldouble import mp_eig, my_eig


data = sio.loadmat("KM.mat")
K = np.array(data['K'], dtype=np.complex256)
M = np.array(data['M'], dtype=np.complex256)

gamma, vecs = sla.eig(K,M)
res = K @ vecs - M @ vecs @ np.diag(gamma.reshape((-1)))
print("Initial residual:", sla.norm(res, 'fro'))

KMI = sla.solve(M.copy(),K.copy())
# gamma, vecs = mp_eig(KMI.copy())
# res = K @ vecs - M @ vecs @ np.diag(gamma)
# print("Improved residual:", sla.norm(res, 'fro'))

gamma, vecs = my_eig(KMI.copy())
res = K @ vecs - M @ vecs @ np.diag(gamma)
print("Improved residual:", sla.norm(res, 'fro'))


# for j in range(vecs.shape[1]):
#     v, sigma = improve(K, M, vecs[:, j], gamma[j])
#     vecs[:, j] = v
#     gamma[j] = sigma










gamma = gamma*1j/1.0
plt.figure()
plt.plot(np.real(gamma), np.imag(gamma), 'ok', markersize=3)
plt.xlabel('Re(γ)')
plt.ylabel('Im(γ)')
plt.xlim((-0.2,1.8))
plt.ylim((-1,0.2))
plt.grid(True)
plt.show()
