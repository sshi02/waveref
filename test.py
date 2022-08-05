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

h1 = 1
h2 = 0.5

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

(a_i, a_r, i_min, i_max, K_r) = waveref.reflection(eta[0], eta[1], dl, dt, h)

print("expected coefficent of reflection: {:f}".format(h2 / h1))
print("calculated coefficent of reflection: {:f}".format(K_r))

fig = plt.figure()
plt.subplot(1, 2, 1)
plt.plot(range(0, 500), a_i, "r-")
plt.subplot(1, 2, 2)
plt.plot(range(0, 500), a_r)

plt.savefig('reflection_test1.png')

print("---")
print("test2: reflection() standing wave")
print("test2 parameters")

dl = 10
dt = 2
h = 8

print("water depth: {:f}".format(h))
print("gauge distance: {:f}".format(dl))
print("reading in station data from wall_case")

with open('./wall_case/sta_0001', 'r') as file:
    eta1 = np.loadtxt(file).flatten()
with open('./wall_case/sta_0002', 'r') as file:
    eta2 = np.loadtxt(file).flatten()
print("---")

(a_i, a_r, i_min, i_max, K_r) = waveref.reflection(eta1, eta2, dl, dt, h)

print("calculated coefficent of reflection: {:f}".format(K_r))
