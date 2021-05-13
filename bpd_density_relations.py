#!/usr/bin/env python

__author__ = "Jan Niederau"
__license__ = "Apache-2.0 License"


import numpy as np


def scaled_temperature(temperature: float=200., pressure: float=981., XNaCl: float=0.5):
    """Volumetric scaling of temperature. After Driesner 2007 II, eq.8 to e.16

    Args:
        temperature (float): scaled temperature for volumetric correlation in degree Celsius. Defaults to 200 °C.
        pressure (float): pressure in bar. Defaults to 470 bar. 4 km depth with a density of 1200 kg/m**3
        XNaCl (float): Mole fraction of NaCl

    Returns:
        [type]: [description]
    """
    # parameters for n1 and n2
    n11 = -54.2958 - 45.7623 * np.exp(-9.44785e-4 * pressure)
    n21 = -2.6142 - 0.000239092 * pressure
    n22 = 0.0356828 + 4.37235e-6 * pressure + 2.0566e-9 * pressure**2
    
    # parameters for n30 and n31
    n300 = 7.60664e6 / (pressure + 472.051)**2
    n301 = -50 - 86.1446 * np.exp(-6.21128e-4 * pressure)
    n302 = 294.318 * np.exp(-5.66735e-3 * pressure)
    n310 = -0.073276 * np.exp(-2.3772e-3 * pressure)
    n311 = -47.2747 + 24.3653 * np.exp(-1.25533e-3 * pressure)
    n312 = -0.278529 - 0.00081381 * pressure

    n1 = n11 * (1- XNaCl) # + n12 * (1 - XNaCl)**2
    n2 = n21 * np.sqrt(XNaCl + n22) # + n23 * XNaCl

    n30 = n300 * (np.exp(n301 * XNaCl) - 1) + n302 * XNaCl
    n31 = n310 * (np.exp(n311 * XNaCl) + n312 * XNaCl)
    Dt = n30 * np.exp(n31*temperature)

    TV = n1 + n2 * temperature + Dt

    return TV