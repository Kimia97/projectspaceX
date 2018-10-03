import numpy as np


def get_air_drag(h, Cd, A, v):
    if h < 0:
        return 0

    return 1/2 * Cd * get_air_resistance(h) * A * v**2


def get_air_resistance(h):
        return get_air_pressure(h) / get_temperature(h) * 3.4855


# Returns the air pressure in Pascal based on height above sea level
def get_air_pressure(h):

    if 0 <= h <= 11000:  # Troposphere
        return 101.29e3 * (get_temperature(h)/288.08)**5.256

    elif 11000 < h <= 25000:  # Lower stratosphere
        return 127.76e3 * np.e**(-0.000157 * h)

    elif 25000 < h:  # Higher stratosphere
        return  2.488e3 * (get_temperature(h)/216.6)**(-11.388)


# Returns the air temperature in Kelvin based on height above sea level
def get_temperature(h):

    if 0 <= h <= 11000:  # Troposphere
        return 288.19 - 0.00649 * h

    elif 11000 < h <= 25000:  # Lower stratosphere
        return 216.69

    elif 25000 < h:  # Higher stratosphere
        return 141.94 + 0.00299 * h

