import numpy as np

def xfer_2_1(f):
    s = 1j * 2*np.pi*f
    R, C = 56e3, 4.7e-9
    omega0 = 1/R/C
    T = omega0/(s+omega0)
    return T

def xfer_2_2(f):
    s = 1j * 2*np.pi*f
    R, C = 56e3, 4.7e-9
    omega0 = 1/R/C
    T = s / (s+omega0)
    return T

def xfer_2_3(f):
    s = 1j * 2*np.pi*f
    R1, R2, C = 220e3, 18e3, 4.7e-9
    RP = R1*R2/(R1+R2)
    omega1 = 1/R1/C
    omegaP = 1/RP/C
    T = (s+omega1) / (s+omegaP)
    return T

def xfer_2_4(f):
    s = 1j * 2*np.pi*f
    R1, R2, C = 18e3, 220e3, 4.7e-9
    omega1 = 1/R1/C
    omega2 = 1/(R1+R2)/C
    T = omega2/omega1 * (s+omega1) / (s+omega2)
    return T

def xfer_T_H(f, R2):
    s = 1j * 2*np.pi*f
    R, C, R1 = 56e3, 4.7e-9, 47e3
    omega0 = 1/R/C
    Q = R2/R1
    T = - s*s / (s*s + (omega0/Q)*s + omega0*omega0)
    return T

def xfer_T_B(f, R2):
    s = 1j * 2*np.pi*f
    R, C, R1 = 56e3, 4.7e-9, 47e3
    omega0 = 1/R/C
    Q = R2/R1
    T = omega0*s / (s*s + (omega0/Q)*s + omega0*omega0)
    return T

def xfer_T_L(f, R2):
    s = 1j * 2*np.pi*f
    R, C, R1 = 56e3, 4.7e-9, 47e3
    omega0 = 1/R/C
    Q = R2/R1
    T = - omega0*omega0 / (s*s + (omega0/Q)*s + omega0*omega0)
    return T