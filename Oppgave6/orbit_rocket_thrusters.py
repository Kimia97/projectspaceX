import numpy as np


def get_rocket_mass(t):
    m1 = 2970000  # kg
    m2 = 680000   # kg
    m3 = 183800   # kg

    dtm1 = (90000/7)   # kg/s
    dtm2 = (22805/18)  # kg/s
    dtm3 = 219         # kg/s

    if t <= 168:
        return m1 - t * dtm1
    elif 160 < t <= 360:
        return m2 - t * dtm2
    elif 360 < t <= 165 + 335:
        return m3 - t * dtm3
    else:
        return m3 - (165 + 335) * dtm3


# Angle from a stick going straight trough the top of the earth
def get_rocket_thrust(t, angle):
    thrust1 = 35100 * 10e3  # N
    thrust2 = 5141 * 10e3   # N
    thrust3 = 1000 * 10e3   # N

    if t <= 168:
        Fy = thrust1 * np.cos(np.deg2rad(angle))
        Fx = thrust1 * np.sin(np.deg2rad(angle))
        return [Fx, Fy]
    elif 160 < t <= 360:
        Fy = thrust2 * np.cos(np.deg2rad(angle))
        Fx = thrust2 * np.sin(np.deg2rad(angle))
        return [Fx, Fy]
    elif 360 < t <= 165 + 335:
        Fy = thrust3 * np.cos(np.deg2rad(angle))
        Fx = thrust3 * np.sin(np.deg2rad(angle))
        return [Fx, Fy]
    else:
        return [0, 0]
