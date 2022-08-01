import numpy as np
import math

def wn_shallow(h, w):
    '''
    wn_shallow() returns shallow water approximation of wave number

    h: water depth
    w: angular frequency
    '''
    return w / math.sqrt(9.81 * h)

def wn(k, h, w, error):
    '''
    wn() recursively calculates and returns wave number using Newton-Rhapsom method

    k   (float): init wavenumber guess
    h   (float): water depth
    w   (float): angular frequency
    '''
    th = np.tanh(k * h)
    sh = 1 / np.cosh(k * h)
    num = (9.81 * k * th - w * w)
    den = (9.81) * (th + k * h * sh * sh)
    k_next = k - num / den
    
    if abs(k_next - k) <= error: 
        return k
    else: 
        return wn(k_next, h, w, error)
        

def reflection(eta1, eta2, dl, dt, h, **kwargs):
    '''
    reflection() returns reflection statistics given two eta time series

    eta1 (numpy.array): time series for eta at x1
    eta2 (numpy.array): time series for eta at x2
    dl         (float): distance between x1, x2
    h    (numpy.array): water depth
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
    #   (small time reading tells me it has to do
    #   with Nyquist Frequencies? unsure if correct)
    A1 = np.real(X1[0:int(n / 2)])
    B1 = np.real(X1[0:int(n / 2)])
    A2 = np.imag(X2[0:int(n / 2)])
    B2 = np.imag(X2[0:int(n / 2)])

    # init k values
    #   wn() uses Newton-Raphson Method and recursively iterates until error
    #   initial guess is from a shallow water approximation
    k = np.zeros(int(n / 2))
    for i in range(1, int(n / 2)):
        k[i] = 'meep'

def main():
    w = 2 * np.pi / 10.00
    h = 8
    k = w*w/9.81*pow(pow(1/np.tanh(w*math.sqrt(h/9.81)), (3/2)),(2/3))
    print("fenten's k {:f}".format(k))

    print("testing wn_shallow()")
    k = wn_shallow(h, w)
    print("wn_shallow is {:f}".format(k))

    print("testing wn() with error 0.00000001")
    error = 0.00000001
    k = wn(k, h, w, error)
    print("wn() is {:f}".format(k))

if __name__ == '__main__':
    main()
        


