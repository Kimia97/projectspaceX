from numpy import sqrt
import time

import numpy as np
import scipy.integrate as integrate

import matplotlib.pyplot as plot
import matplotlib.animation as animation

class Orbit:
    """
    
    Orbit Class

    init_state is [t0,x0,vx0,y0,vx0],
    where (x0,y0) is the initial position
    , (vx0,vy0) is the initial velocity
    and t0 is the initial time
    """
    def __init__(self,
                 init_state=[0, 0, 1, 2, 0],
                 G=3,
                 m_earth=5,   # Jordens masse er 5.9736 * 10^24 Kg
                 m_moon=5):  # Månens masse er 7.3477 * 10^22 Kg
        self.GravConst = G
        self.m_earth = m_earth
        self.m_moon = m_moon
        self.state = np.asarray(init_state, dtype='float')
    
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
        U = -G*m1*m2/sqrt(x**2+y**2)
        K = m1*(vx**2+vy**2)/2
        return K+U

    def time_elapsed(self):
        return self.state[0]

    def get_velocity(self):
        return [self.state[2], self.state[4]]

    def get_position(self):
        return [self.state[1], self.state[3]]


    def step(self, h):
        """Uses the trapes method to calculate the new state after h seconds."""
        x = self.state
        s1 = self.ydot(x)
        s2 = self.ydot(x+h*s1)
        self.state = x+h*(s1+s2)/2
    
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
        dist = sqrt((px2-px1)**2+(py2-py1)**2)
        z = np.zeros(5)
        z[0] = 1
        z[1] = vx1
        z[2] = (Gm2*(px2-px1))/(dist**3)
        z[3] = vy1
        z[4] = (Gm2*(py2-py1))/(dist**3)

        return z


# make an Orbit instance

orbit = Orbit([0.0, 0.0, 1.2, 10, 0.0])
dt = 1./30  # 30 frames per second

# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-12, 12), ylim=(-12, 12))

line1, = axes.plot([], [], 'o-r', linewidth=10, markersize=4)  # Dette er månen
line2, = axes.plot([], [], 'o-b', linewidth=10, markersize=16)  # Dette er jorden

time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
energy_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)
velocity_text = axes.text(0.02, 0.85, '', transform=axes.transAxes)
position_text = axes.text(0.02, 0.80, '', transform=axes.transAxes)



def init():
    """initialize animation"""
    line1.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    energy_text.set_text('')
    return line1, line2, time_text, energy_text, velocity_text, position_text


def animate(i):
    """perform animation step"""
    global orbit, dt
    orbit.step(dt)
    line1.set_data(*orbit.position())
    line2.set_data([0.0, 0.0])
    time_text.set_text('time = %.1f' % orbit.time_elapsed())
    energy_text.set_text('energy = %.3f J' % orbit.energy())
    velocity_text.set_text('velocity = %.3f x' % orbit.get_velocity()[0]+ ', %.3f y' % orbit.get_velocity()[1])
    position_text.set_text('x = %.3f ' % orbit.get_position()[0]+ ', %.3f y' % orbit.get_position()[1])

    return line1, line2, time_text, energy_text, velocity_text, position_text


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 1000 * dt - (t1 - t0)

anim = animation.FuncAnimation(fig,             # figure to plot in
                               animate,         # function that is called on each frame
                               frames=6000,      # total number of frames
                               interval=delay,  # time to wait between each frame.
                               repeat=False,
                               blit=True,
                               init_func=init   # initialization
                               )

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
# anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plot.show()

print('DONE')
