from numpy import sqrt
from Oppgave2 import oppgave2 as rk45
import time
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.animation as animation


class Orbit:
    """
    Orbit Class

    init_state is [t0,x0,vx0,y0,vy0],
    where (x0,y0) is the initial position
    , (vx0,vy0) is the initial velocity
    and t0 is the initial time
    """

    def __init__(self,
                 init_state=[0, 0, 1, 2, 0],
                 G=6.67e-11,
                 m_earth = 5.97e24):        # m_earth=6e24):
        self.GravConst = G
        self.m_earth = m_earth  # Jordens masse er 5.9736 * 10^24 Kg
        self.m_moon = 7.3477e22  # Månens masse er 7.3477 * 10^22 Kg, 0.012 x Jordens
        self.state = np.asarray(init_state, dtype='float')
        self.orbit_time = 0.0
        self.xy = [[0], [0.4055e9]]

    def position(self):
        """compute the current x,y positions of the pendulum arms"""
        x = self.state[1]
        y = self.state[3]
        return x, y

    def energy(self):
        x = self.state[1]
        y = self.state[3]
        vx = self.state[2]
        vy = self.state[4]
        m1 = self.m_earth
        m2 = self.m_moon
        G = self.GravConst
        U = -G * m1 * m2 / sqrt(x ** 2 + y ** 2)
        K = m1 * (vx ** 2 + vy ** 2) / 2
        return K + U

    def time_elapsed(self):
        return self.state[0]

    def get_velocity(self):
        return [self.state[2], self.state[4]]

    def get_position(self):
        return [self.state[1], self.state[3]]

    def get_orbit_time(self):
        if (self.state[4] - 4) <= 0 <= self.state[4]:
            self.orbit_time = self.time_elapsed() / 86400

        return self.orbit_time

    def step(self, h):
        x = self.state
        s1 = self.ydot(x)
        s2 = self.ydot(x + h * s1)
        self.state = x + h * (s1 + s2) / 2

    def ydot(self, x):
        G = self.GravConst
        m2 = self.m_earth
        Gm2 = G * m2

        px2 = 0
        py2 = 0
        px1 = x[1]
        py1 = x[3]
        vx1 = x[2]
        vy1 = x[4]
        dist = sqrt((px2 - px1) ** 2 + (py2 - py1) ** 2)
        z = np.zeros(5)

        z[0] = 1

        z[1] = vx1
        z[2] = (Gm2 * (px2 - px1)) / (dist ** 3)
        z[3] = vy1
        z[4] = (Gm2 * (py2 - py1)) / (dist ** 3)

        self.xy[0].append(self.get_position()[0])
        self.xy[1].append(self.get_position()[1])

        return z


# make an Orbit instance
orbit = Orbit([0.0, 0.0, 970, 0.4055e9, 0.0])

dt = 1.0/30.0  # 30 frames per second

# Use runge-kutta 4/5
rk45 = rk45.RungeKuttaFehlberg54(orbit.ydot, 5, dt, 05e-14)

# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-4.5e8, 4.5e8), ylim=(-4.5e8, 4.5e8))

line1, = axes.plot([], [], 'ro', linewidth=1, markersize=10)  # Dette er månen
line1_2, = axes.plot([], [], 'r--', linewidth=0.5)  # Dette er linjen som viser hvor månen har vært
line2, = axes.plot([], [], 'bo', linewidth=10, markersize=30)  # Dette er jorden

time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
velocity_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)
position_text = axes.text(0.02, 0.85, '', transform=axes.transAxes)
orbit_time_text = axes.text(0.02, 0.80, '', transform=axes.transAxes)


def init():
    """initialize animation"""

    line1.set_data([], [])
    line1_2.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    # energy_text.set_text('')
    return line1, line1_2, line2, time_text, velocity_text, position_text, orbit_time_text


def animate(i):
    """perform animation step"""
    global orbit, dt
    orbit.state, E = rk45.safeStep(orbit.state)
    line1.set_data(*orbit.position())
    line1_2.set_data(orbit.xy)
    line2.set_data([0.0, 0.0])
    time_text.set_text('time = %.3f dager' % (orbit.time_elapsed() / 86400))
    # energy_text.set_text('energy = %.3f J' % orbit.energy())
    velocity_text.set_text('velocity = %.3f x km/s' % (orbit.get_velocity()[0]/1e3) + ', %.3f y km/s' % (orbit.get_velocity()[1]/1e3))
    position_text.set_text('x = %.3f*10e6 km' % (orbit.get_position()[0]/1e9) + ', y = %.3f*10e6 km' % (orbit.get_position()[1]/1e9))
    orbit_time_text.set_text('orbit time  = %.3f dager' % orbit.get_orbit_time())

    return line1, line1_2, line2, time_text, velocity_text, position_text, orbit_time_text


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 1000 * dt - (t1 - t0)

anim = animation.FuncAnimation(fig,  # figure to plot in
                               animate,  # function that is called on each frame
                               frames=1000,  # total number of frames
                               interval=delay,  # time to wait between each frame.
                               repeat=False,
                               blit=True,
                               init_func=init  # initialization
                               )

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
# anim.save('orbit.mp4', fps=35, extra_args=['-vcodec', 'libx264'])

plot.show()

print('DONE')
