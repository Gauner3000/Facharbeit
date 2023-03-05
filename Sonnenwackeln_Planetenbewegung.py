import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter







class Planet:

    def __init__(self, x, y, x_vel, y_vel, mass, color):
        self.x = x
        self.y = y
        self.x_vel = x_vel
        self.y_vel = y_vel
        self.mass = mass
        self.color = color

        self.orbit = []

#sun = Planet(x-Koordiante (vom Ursprung), y-Koordinate (vom Ursprung), x-Geschwindigkeit (in m/s), y-Geschwindigkeit (in m/s), Masse (kg), Farbe)
sun = Planet(0, 0, 0, 0, 1.989e30, "yellow") #Sonne
earth = Planet(149.6e9, 0, 0, 30e8, 5.972e24, "blue") #Erde
Objekt = Planet(227.9e9, 0, 0, 24e8, 6.39e23, "orange") #Mars





fig = plt.figure(figsize=(6.4*1.5, 4.8*1.5) )
fig.tight_layout(pad=4.0)
#---------------------------------------------------------------
ax = fig.add_subplot(2, 2, 1)
ax.set(xlim=(-4.5e11, 4e11), ylim=(-4e11, 4e11)) #Bewegung Planeten
ax.set_title("Bewegung der Körper")
ax.set_xlabel("m")
ax.set_ylabel("m")

ln1, = plt.plot([], [], color=sun.color, marker="o")
ln2, = plt.plot([], [], color=earth.color, marker="o")
ln3, = plt.plot([], [], color=Objekt.color, marker="o")

#------------------------------------------------------------
bx = fig.add_subplot(2, 2, 2)
bx.set(xlim=(-0.5e8, 0.5e8), ylim=(-0.5e8, 0.5e8)) #Bewegung Sonne
bx.set_title("Bewegung Sonne")
bx.set_xlabel("m")
bx.set_ylabel("m")

lb1, = plt.plot([], [], color=sun.color, marker="o")





def dSdt(t, S):

    x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3 = S
    r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    r13 = np.sqrt((x3-x1)**2 + (y3-y1)**2)
    r23 = np.sqrt((x2-x3)**2 + (y2-y3)**2)


    return [ vx1,
            vy1,
            vx2,
            vy2,
            vx3,
            vy3,
            earth.mass/r12**3 * (x2-x1) + Objekt.mass/r13**3 * (x3-x1), #mass 1
            earth.mass/r12**3 * (y2-y1) + Objekt.mass/r13**3 * (y3-y1),
            sun.mass/r12**3 * (x1-x2) + Objekt.mass/r23**3 * (x3-x2), #mass 2
            sun.mass/r12**3 * (y1-y2) + Objekt.mass/r23**3 * (y3-y2),
            sun.mass/r13**3 * (x1-x3) + earth.mass/r23**3 * (x2-x3), #mass 3
            sun.mass/r13**3 * (y1-y3) + earth.mass/r23**3 * (y2-y3)
           ]


t = np.linspace(0, 2000, 1000)

sol = solve_ivp(dSdt, (0, 2000), y0=[sun.x, sun.y, earth.x, earth.y, Objekt.x, Objekt.y,
                       sun.x_vel, sun.y_vel, earth.x_vel, earth.y_vel, Objekt.x_vel, Objekt.y_vel],
                method = 'DOP853', t_eval=t, rtol=1e-10, atol=1e-13)


t = sol.t
x1 = sol.y[0]
y1 = sol.y[1]
x2 = sol.y[2]
y2 = sol.y[3]
x3 = sol.y[4]
y3 = sol.y[5]
#plt.plot(t, x1)
#plt.show()


def animate(i):
    ln1.set_data(x1[i], y1[i])
    ln2.set_data(x2[i], y2[i])
    ln3.set_data(x3[i], y3[i])

    lb1.set_data(x1[i], y1[i])

    return  ln1, ln2, ln3, lb1



ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
plt.show()