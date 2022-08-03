import numpy as np
from matplotlib import pyplot as plt
import waveref

print("==============================================")
print("testing wave number functions")
print("==============================================")
print("test parameters")
w = 2 * np.pi / 10.00
h = 8
k = 2 * np.pi / 83.81
print("omega: {:f}".format(w))
print("water depth: {:f}".format(h))
print("wave number: {:f}".format(k))
print("---")

print("testing wn_shallow()")
k = waveref.wn_shallow(h, w)
print("wn_shallow() is {:f}".format(k))
print("---")

error = 0.000001
print("testing wn() with error {:f}".format(error))
k = waveref.wn(k, h, w, error)
print("wn() is {:f}".format(k))
print("---")

print("==============================================")
print("testing reflection()")
print("==============================================")
print("test1: reflection() opposing cosine waves test")
print("test1 parameters")
w = 2 * np.pi / 10.00
h = 8
k = 2 * np.pi / 83.81
T = np.transpose(np.array(range(1000)))
dt = 1
x = np.array([100, 110])
dl = x[1] - x[0]

h1 = 0.5
h2 = 1

kx = k * x
eta1 = np.zeros((2, 1000))
eta2 = np.zeros((2, 1000))
eta1[0] = h1 / 2 * np.cos(kx[0] - T * w)
eta1[1] = h1 / 2 * np.cos(kx[1] - T * w)
eta2[0] = h2 / 2 * np.cos(kx[0] + T * w)
eta2[1] = h2 / 2 * np.cos(kx[1] + T * w)

eta = eta1 + eta2

print("omega: {:f}".format(w))
print("water depth: {:f}".format(h))
print("wave number: {:f}".format(k))
print("gauge distance: {:f}".format(dl))
print("incident wave height: {:f}".format(h1))
print("reflected wave height: {:f}".format(h2))
print("---")

(a_i, a_r) = waveref.reflection(eta[0], eta[1], dl, dt, h)

fig = plt.figure()
plt.subplot(1, 2, 1)
plt.plot(range(0, 500), a_i, "r-")
plt.subplot(1, 2, 2)
plt.plot(range(0, 500), a_r)

plt.savefig('reflection_test1.png')