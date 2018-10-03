from numpy import sqrt
from Oppgave2 import oppgave2 as rk45
from Oppgave6 import orbit_rocket_thrusters as ort
from Oppgave4 import rocket
from Oppgave5 import atmosphere
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
                 m_earth=5.9736e24):
        self.GravConst = G
        self.m_earth = m_earth  # Jordens masse er 5.9736 * 10^24 Kg
        self.m_rocket = rocket.get_rocket_mass(0)
        self.state = np.asarray(init_state, dtype='float')
        self.xy = [[0], [6371e3]]
        self.last_time = -1
        self.acceleration = 0
        self.rocket_angle = 0

    def position(self):
        """compute the current x,y positions of the pendulum arms"""
        x = self.state[1]
        y = self.state[3]
        return x, y

    def distance_above_sea_level(self):
        x = self.state[1]
        y = self.state[3]
        earth_radius = 6371e3
        distance_from_earth_center = sqrt(x*x + y*y)

        return distance_from_earth_center - earth_radius

    def get_angle(self):
        return 90 + np.rad2deg(np.arctan(self.state[1]/self.state[3]))

    def time_elapsed(self):
        return self.state[0]

    def get_abs_vel(self):
        return np.sqrt(self.state[2]**2 + self.state[4]**2)

    def get_velocity(self):
        return [self.state[2], self.state[4]]

    def get_acceleration(self):
        return self.acceleration

    def get_position(self):
        return [self.state[1], self.state[3]]

    def step(self, h):
        """Uses the trapezoid method to calculate the new state after h seconds."""
        x = self.state
        s1 = self.ydot(x)
        s2 = self.ydot(x + h * s1)
        self.state = x + h * (s1 + s2) / 2

    def ydot(self, x):
        G = self.GravConst
        m2 = self.m_earth
        Gm2 = G * m2
        t = x[0]
        px2 = 0
        py2 = 0
        px1 = x[1]
        py1 = x[3]
        vx1 = x[2]
        vy1 = x[4]
        dist = sqrt((px2 - px1) ** 2 + (py2 - py1) ** 2)
        z = np.zeros(5)

        # Force from gravity on rocket divided by rocket mass
        Fg_x = (Gm2 * (px2 - px1)) / (dist ** 3)
        Fg_y = (Gm2 * (py2 - py1)) / (dist ** 3)

        # Area of the rocket used in the drag calculations
        a1 = np.pi * 5.55**2

        rocket_mass = ort.get_rocket_mass(t)

        # Used to calculate the angle of the speed relative to starting angle
        if vy1 == 0 or vx1 == 0:
            speed_angle = 0
        elif vy1 < 0 and vx1 >= 0:
            speed_angle = np.pi + np.arctan(vx1 / vy1)
        elif vy1 < 0 and vx1 <=0:
            speed_angle = np.pi/2 + np.arctan(vx1 / vy1)
        elif vy1 > 0 and vx1 <= 0:
            speed_angle = np.pi + np.arctan(vx1 / vy1)
        elif vy1 > 0 and vx1 >= 0:
            speed_angle = 0 + np.arctan(vx1 / vy1)

        # After 130 seconds increase angle until speed is 90 degrees from starting angle
        if t > 130:
            if np.rad2deg(speed_angle) >= 90:
                self.rocket_angle = self.get_angle()
            elif round(t) != self.last_time and round(t) > self.last_time:
                self.rocket_angle += 3
                self.last_time = round(t)


        # Force from thrusters on rocket divided by rocket
        F = ort.get_rocket_thrust(t, self.rocket_angle)
        Fx, Fy = F[0], F[1]

        # Switch off rocket thrusters after 175 seconds
        if t >= 175:
            rocket_mass = ort.get_rocket_mass(175)
            Fx, Fy = 0, 0


        # Force from air drag on rocket divided by rocket mass
        abs_vel = np.sqrt(vx1**2 + vy1**2)
        Fd = atmosphere.get_air_drag(self.distance_above_sea_level(), 0.2, a1, abs_vel)/rocket_mass
        Fd_x = np.sin(speed_angle) * Fd
        Fd_y = np.cos(speed_angle) * Fd


        self.acceleration = Fy/rocket_mass - (Fg_y + Fd_y)

        z[0] = 1
        z[1] = vx1
        z[2] = Fx/rocket_mass + (Fg_x - Fd_x)
        z[3] = vy1
        z[4] = Fy/rocket_mass + (Fg_y - Fd_y)

        self.xy[0].append(self.get_position()[0])
        self.xy[1].append(self.get_position()[1])

        return z


# make an Orbit instance
orbit = Orbit([0.0, 0.0, 0, 6371e3, 0.0])

plotScale = 10e6  # meters

dt = 60

# Use runge-kutta 4/5
rk45 = rk45.RungeKuttaFehlberg54(orbit.ydot, 5, dt, 05e-14)

# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-plotScale, plotScale), ylim=(-plotScale, plotScale * 1.5))

earth = plot.Circle((0, 0), 6371e3, color='blue', alpha=0.2)
axes.add_artist(earth)

line1, = axes.plot([], [], 'ro', linewidth=1, markersize=1)  # Dette er raketten
line1_2, = axes.plot([], [], 'r--', linewidth=0.5)  # Dette er linjen som viser hvor månen har vært
# line2, = axes.plot([], [], 'bo', linewidth=10, markersize=10)  # Dette er jorden

time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
velocity_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)
position_text = axes.text(0.02, 0.85, '', transform=axes.transAxes)
height_text = axes.text(0.02, 0.80, '', transform=axes.transAxes)


def init():
    """initialize animation"""
    line1.set_data([], [])
    line1_2.set_data([], [])
    #line2.set_data([], [])
    time_text.set_text('')
    # energy_text.set_text('')
    return line1, line1_2, time_text, velocity_text, position_text, height_text


def animate(i):
    """perform animation step"""
    global orbit, dt, time_0, time_difference

    orbit.state, E = rk45.safeStep(orbit.state)
    line1.set_data(*orbit.position())
    line1_2.set_data(orbit.xy)
    time_text.set_text('time = %.f sek' % orbit.time_elapsed())
    velocity_text.set_text('velocity = %.3f x km/s' % (orbit.get_abs_vel()/1e3))
    height_text.set_text('moh. = %.3f km' % (orbit.distance_above_sea_level() / 1e3))
    position_text.set_text('x = %.3f km' % (orbit.get_position()[0]) + ', y = %.3f km' % ((orbit.get_position()[1] - 6371e3)/1e3))
    return line1, line1_2, time_text, velocity_text, position_text, height_text


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = (1000 * dt - (t1-t0)) * 0

anim = animation.FuncAnimation(fig,  # figure to plot in
                               animate,  # function that is called on each frame
                               frames=300000,  # total number of frames
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
# anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plot.show()

print('DONE')
