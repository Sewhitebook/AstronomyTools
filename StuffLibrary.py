'''This file contains a bunch of different stuff that I have already made, so I may simply pull it from here without having to remake it every time'''

import numpy as np
from scipy import constants

"""Calculus"""

#Numerical Differentiation
def deriv(x, y, xrange): #x, y are the data used to find the derivatives. xrange is the points at which the derivative is desired.
    slopes = []
    for i in range(1, len(y)-1):
        central_diff = (y[i-1] - y[i+1])/(x[i-1] - x[i+1])
        slopes.append(central_diff)

    interp_deriv = interpolate.interp1d(x[1:-1], slopes, kind="cubic", bounds_error=False, fill_value="extrapolate")
    derivatives = interp_deriv(xrange)
    return derivatives

#Secant Method of Differentiating a known function to any desired precision.
def Sec_Method(func, x_guess, tol, maxiter, h, *func_args):
    x = x_guess
    # tol = 1e-7
    # maxiter = 50
    # h = 1e-4
    for i in range(maxiter):
        deriv = (func(x + h, *func_args) - func(x - h, *func_args))/(2*h)
        increment = (0.5 - func(x, *func_args))/deriv
        x += increment
        if np.abs(increment) < tol:
            break
    return x

#Trapezoidal Integration
def Trap(f, x1, x2, N):
    h = (x2 - x1)/N
    t = np.linspace(x1, x2, N + 1)
    weights = np.ones(len(t))*h
    weights[0] = weights[-1] = h/2
    def ret(*args):
        return np.sum(weights*f(args))
    return ret

"""Distributions and Statistics"""

#Gaussian Normal
def Gaussian(x, mean, std):
    pi = np.pi
    a = 1/(std*np.sqrt(2*pi))
    b = np.exp(-1/2*((x-mean)/std)**2)
    ret = a*b
    return ret

#Maxwell-Boltzmann
def Boltzmann(frequency, temperature, mass):
    k = constants.k
    x = (mass/(2*np.pi*k*temperature))**(3/2)
    y = np.exp(-mass*frequency**2/(2*k*temperature))
    ret = 4*np.pi*(frequency**2)*x*y
    return ret

#General Covariance Functions

#constant
def Constant(x1, x2, a):
    _x1, _x2 = np.meshgrid(x1,x2)
    r = a**2*(_x1 - _x2)
    return r

#Squared Exoponential
def Squared_Exp(x1, x2, a, h):
    _x1, _x2 = np.meshgrid(x1, x2)
    r = a**2*np.exp(-(_x1 - _x2)**2/(2*h**2))
    return r

#Ornstein–Uhlenbeck
def Orn_Uhl(x1, x2, a, h):
    _x1, _x2 = np.meshgrid(x1, x2)
    r = a**2*np.exp(-(np.abs(_x1 - _x2)/h))
    return r

#Periodic
def Periodic(x1, x2, a, h):
    _x1, _x2 = np.meshgrid(x1, x2)
    r = a**2*np.exp(-(2*np.sin(_x1/2 - _x2/2)**2/h**2))
    return r

#Generalized Gaussian Process Regression
def GPR(params, sig_data, x, model_domain, y, cov=Constant, scaling=1):
    a = params[0]
    h = params[1]

    Cov_data = scaling*np.identity(len(x))*sig_data**2
    Kinv = np.linalg.inv(scaling*cov(x, x, a, h) + Cov_data)

    y_mean = np.linalg.multi_dot([scaling*cov(x, model_domain, a, h), Kinv, y])

    yvar = scaling*cov(model_domain, model_domain, a, h) 
    yvar -= np.linalg.multi_dot([scaling*cov(x, model_domain, a, h), Kinv, scaling*cov(model_domain, x, a, h)])
    gp_err = np.sqrt(np.diag(yvar))
    return y_mean, gp_err

"""Astronomy"""

    rad_Earth = 6.371e6
    rad_Jupiter = 6.991e7
    rad_Sun = 6.957e8

    mass_Earth = 5.972e24
    mass_Jupiter = 1.898e27
    mass_Sun = 1.988e30

    #For the converters "To" and "Given" must each take the value of "Meters", "Earth", "Jupiter", or "Sun"
    #Radius Converter
    def Radius_Conv(radius, To="Jupiter", Given="Meters"):
        conv_to = "rad_" + To
        conv_from = "rad_" + Given

        if Given == "Meters":
            r = radius/conv_to

        if Given != "Meters":
            r_met = radius*conv_from
            r = r_met/conv_to
        return r
    
    #Radius Converter
    def Mass_Conv(mass, To="Jupiter", Given="Meters"):
        conv_to = "mass_" + To
        conv_from = "mass_" + Given

        if Given == "Meters":
            m = mass/vars(conv_to)

        if Given != "Meters":
            m_met = mass*conv_from
            m = m_met/conv_to
        return m


    #Orbital Period 
    def Orb_Period(mass_star, sm_axis):
        T = 2*np.pi*np.sqrt(sm_axis**3/(mass_star*6.674e-11))
        return T


