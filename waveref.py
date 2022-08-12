import numpy as np
from math import sqrt, ceil     # better specification -> small optimization
import warnings

def wn_shallow(h, w):
    '''
    wn_shallow() returns shallow water approximation of wave number

    h: water depth
    w: angular frequency
    '''
    return w / sqrt(9.81 * h)

def wn(k, h, w, error):
    '''
    wn() recursively calculates and returns wave number 
    using Newton-Rhapsom method

    parameters
    k       (float): init wavenumber guess
    h       (float): water depth
    w       (float): angular frequency
    error   (float): margin of between guesses
    '''
    warnings.filterwarnings("ignore")       # cosh will reach some diverging condition, this is ok
    th = np.tanh(k * h)                     # encaps tanh operations
    sh = 1 / np.cosh(k * h)                 # encaps sech operations
    num = (9.81 * k * th - w ** 2)          # numerator
    den = (9.81) * (th + k * h * sh ** 2)   # denominator 
    k_next = k - num / den
    
    if abs(k_next - k) <= error:    # base case:
        return k
    else:                           # recursive case:
        return wn(k_next, h, w, error)
        

def reflection(eta1, eta2, dl, dt, h, **kwargs):
    '''
    reflection() returns reflection statistics given two eta time series

    parameters:
    eta1    (numpy.array): time series for eta at x1
    eta2    (numpy.array): time series for eta at x2
    dl      (float): distance between x1, x2
    h       (float): water depth
    **kwargs:
    f_min =     (int): minimum resolvable frequency
    f_max =     (int): maximum resolvable frequency


    return
    a_i     (numpy.array): amplitude of incident wave as frequency series
    a_r     (numpy.array): amplitude of reflected wave as frequency series
    i_min   (int): lower array bound given bounded frequency
    i_max   (int): upper array bound given bounded frequency
    K_r     (float): coefficient of reflection
    '''

    n = len(eta1)       # length of time series
    
    # input param checks
    if not len(eta2) == n:
        print("error: eta1 eta2 length not equal")
        return
    
    # even t-series check, because the FFT of eta
    #   will be normalized and indexed bound by N/2
    if not n % 2 == 0:
        print("warning: time series truncated to even length")
        n -= 1

    # FFT of eta1, eta2, normalized relative to N/2
    (X1, X2) = (np.fft.fft(eta1, n) / (n/2), 
                np.fft.fft(eta2, n) / (n/2))

    # init A1, A2, B1, B2 by collecting sin/cos and truncating
    #   half time-series due to symmetry 
    A1 = np.real(X1[0:int(n / 2)])
    B1 = np.imag(X1[0:int(n / 2)])
    A2 = np.real(X2[0:int(n / 2)])
    B2 = np.imag(X2[0:int(n / 2)])

    # from matplotlib import pyplot as plt
    # plt.figure()
    # plt.plot(range(0, n), X1, "r-")
    # plt.savefig('nontrucated.png')


    # init k values
    #   wn() uses Newton-Raphson Method and     recursively iterates 
    #       until desired closeness
    #   initial guess is from a shallow water approximation
    k = np.zeros(int(n / 2))
    w = np.zeros(int(n / 2))
    for i in range(1, int(n / 2)):
        w[i] = 2 * np.pi * i / dt / n
        k[i] = wn(wn_shallow(h, w[i]), h, w[i], 0.00001)
        #k[i] = w*w/9.81*pow(pow(1/np.tanh(w*sqrt(h/9.81)), (3/2)),(2/3)) # fentons approximation

    # amplitude calculation
    # equation (5), Goda 76
    a_i = np.zeros((int(n/2)))
    a_r = np.zeros((int(n/2)))
    for i in range(int(n/2)):
        den = 2 * abs(np.sin(k[i] * dl))        # init denominator calculation
        if den == 0:                            # diverging condition
            a_i[i] = np.NaN
            a_r[i] = np.NaN
        else: 
            sqr1 = A2[i] - A1[i] * np.cos(k[i] * dl) - B1[i] * np.sin(k[i] * dl)
            sqr2 = B2[i] + A1[i] * np.sin(k[i] * dl) - B1[i] * np.cos(k[i] * dl)
            sqr3 = A2[i] - A1[i] * np.cos(k[i] * dl) + B1[i] * np.sin(k[i] * dl) 
            sqr4 = B2[i] - A1[i] * np.sin(k[i] * dl) - B1[i] * np.cos(k[i] * dl) 
            a_r[i] = np.sqrt(np.square(sqr1) + np.square(sqr2)) / den
            a_i[i] = np.sqrt(np.square(sqr3) + np.square(sqr4)) / den  
    
    # frequency limits
    #   equation (7), Goda 76
    if 'f_min' in kwargs:               # lower bound
        f_min = float(kwargs.get('f_min'))
    else: 
        f_min = 0.05 
    if 'f_max' in kwargs:               # upper bound
        f_max = float(kwargs.get('f_max'))
    else:
        f_max = 0.45

    # array limits
    i_min = int(ceil(f_min * n * dt))       # lower bound
    i_max = int(f_max * n * dt)             # upper bound

    # energy calculation, bounded by above limits
    #   equation (8), Goda 76
    e_i = np.sum(np.square(w[i_min:i_max]) * np.square(a_i[i_min:i_max]) / 2)
    e_r = np.sum(np.square(w[i_min:i_max]) * np.square(a_r[i_min:i_max]) / 2)
    
    # coefficent of reflection calculation
    #   equation (9), Goda 76
    K_r = sqrt(e_r / e_i)

    return (a_i, a_r, i_min, i_max, K_r)

def main():
    print("waveref.py main() invoked")

if __name__ == '__main__':
    main()