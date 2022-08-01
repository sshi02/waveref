import numpy as np
import math
import waveref

print("==============================================")
print("testing wave number functions")
print("==============================================")
print("test parameters")
w = 2 * np.pi / 10.00
h = 8
print("omega: {:f}".format(w))
print("water depth: {:f}".format(h))

print("---")
k = w*w/9.81*pow(pow(1/np.tanh(w*math.sqrt(h/9.81)), (3/2)),(2/3))
print("fenten's k is {:f}".format(k))
print("---")

print("testing wn_shallow()")
k = waveref.wn_shallow(h, w)
print("wn_shallow is {:f}".format(k))
print("---")

error = 0.000001
print("testing wn() with error {:f}".format(error))
k = waveref.wn(k, h, w, error)
print("wn() is {:f}".format(k))
print("---")