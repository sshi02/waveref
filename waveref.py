from matplotlib import pyplot
import numpy as np
import math

def wn(k, h, T, error):
    '''
    wn() recursively calculates wavelength using Newton-Rhapsom method

    L   (float): init wavelength guess
    h   (float): water depth
    T   (float): period
    '''


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
        k[i] = 2*np.pi / wn(2 * np.pi / (math.sqrt(9.8 * h) / (n * dt / i)), h, n * dt / i)

def main():
    print("running main() waveref.py")

if __name__ == '__main__':
    main()
        


