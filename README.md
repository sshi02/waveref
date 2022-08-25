# waveref.py
waveref.py is a module that provides methods for evaluating 1-D wave reflection. Largely dependent on Goda and Suzuki 76.

### waveref.**reflection**(*eta1, eta2, dl, dt, h, kwargs*)
*reflection() returns reflection statistics given two eta time series*  
#### parameters
- eta1    (numpy.array): time series for eta at x1
- eta2    (numpy.array): time series for eta at x2
- h       (float): water depth  
- f_min=     (int): minimum resolvable frequency
- f_max=     (int): maximum resolvable frequency
#### return
- a_i     (numpy.array): amplitude of incident wave as frequency series
- a_r     (numpy.array): amplitude of reflected wave as frequency series
- i_min   (int): lower array bound given bounded frequency
- i_max   (int): upper array bound given bounded frequency
- e_i     (float): energy of incident wave
- e_r     (float): energy of reflected wave
- K_r     (float): coefficient of reflection

#### example usage
simple half reflection test [test.py --> test 1] - the following code creates two time series of 'measured' surface elevation at two locations for half reflection. See test1.png for amplitude comparison.
```
import waveref
import numpy as np

# create a fake time series for incident and reflected components
w = 2 * np.pi / 10.00         # angular frequency, f = 10
h = 8                         # water depth, 8 (m)
k = 2 * np.pi / 83.81         # angular wavenumber, wavelength = 83.81
dt = 0.1                      # measure timestep by gauges, 0.1 (s)
T = np.transpose(np.array(range(3600)) * dt)    # time series 
x = np.array([100, 110])      # fake gauge positions
dl = x[1] - x[0]              # distance between gauges
h_i = 1                       # incident wave height
h_r = 0.5                     # reflected wave height
                                                                          # superimposing two opposing waves
eta1 = h1 / 2 * np.cos(kx[0] - T * w) + h2 / 2 * np.cos(kx[0] + T * w)    # at gauge 1 (incident) (x = 100m)
eta2 = h1 / 2 * np.cos(kx[1] - T * w) + h2 / 2 * np.cos(kx[1] + T * w)    # at gauge 2 (incident) (x = 110m)

(a_i, a_r, i_min, i_max, e_i, e_r, K_r) = waveref.reflection(eta1, eta2, dl, dt, h)
```

### waveref.**wn_shallow**(*h, w*)
wn_shallow() returns shallow water approximation of wavenumber
#### parameters
- h     (float): water depth
- w     (float): angular frequency
#### return
- k     (float): wavenumber

### waveref.**wn**(*k, h, w, error, c*)
wn() recursively calculates and returns wave number using Newton-Rhapsom method
#### parameters
- k     (float): wavenumber
- h     (float): water depth
- w     (float): angular frequency
- error (float): margin between guesses
- c     (float): recursive counter (usually set to 0)
#### return
- k     (float): wavenumber
