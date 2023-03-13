import scipy as sci
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import scipy.integrate

#Definitionen
G=6.67408e-11
m_nd=1.989e+30 #Masse der Sonne
r_nd=5.326e+12
v_nd=30000
t_nd=79.91*365*24*3600*0.51
K1=G*t_nd*m_nd/(r_nd**2*v_nd)
K2=v_nd*t_nd/r_nd

#Definition der Massen
m1=1.1 #Alpha Centauri A
m2=0.907 #Alpha Centauri B
m3=1.0 #Dritter Stern
#Definition der Anfangs-Positionen
r1=np.array([-0.5,0,0], dtype="float64")
r2=np.array([0.5,0,0], dtype="float64")
r3=np.array([0,1,0], dtype="float64")

#Definition der Anfangs-Geschwindigkeiten
v1=np.array([0.01,0.01,0], dtype="float")
v2=np.array([-0.05,0,-0.1], dtype="float64")
v3=np.array([0,-0.01,0], dtype="float64")

#Updaten der COM Formeln
r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)
v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)


#Bewegungsgleichungen 
def ThreeBodyEquations(w,t,G,m1,m2,m3):
    r1=w[:3]
    r2=w[3:6]
    r3=w[6:9]
    v1=w[9:12]
    v2=w[12:15]
    v3=w[15:18]
    r12=sci.linalg.norm(r2-r1)
    r13=sci.linalg.norm(r3-r1)
    r23=sci.linalg.norm(r3-r2)
    
    dv1bydt=K1*m2*(r2-r1)/r12**3+K1*m3*(r3-r1)/r13**3
    dv2bydt=K1*m1*(r1-r2)/r12**3+K1*m3*(r3-r2)/r23**3
    dv3bydt=K1*m1*(r1-r3)/r13**3+K1*m2*(r2-r3)/r23**3
    dr1bydt=K2*v1
    dr2bydt=K2*v2
    dr3bydt=K2*v3
    r12_derivs=np.concatenate((dr1bydt,dr2bydt))
    r_derivs=np.concatenate((r12_derivs,dr3bydt))
    v12_derivs=np.concatenate((dv1bydt,dv2bydt))
    v_derivs=np.concatenate((v12_derivs,dv3bydt))
    derivs=np.concatenate((r_derivs,v_derivs))
    return derivs


init_params=np.array([r1,r2,r3,v1,v2,v3])
init_params=init_params.flatten() #Erstellen eines 1D Array
time_span=np.linspace(0,20,500) #20 Perioden und 500 Punkte

#Integrieren der Funktion
three_body_sol=sci.integrate.odeint(ThreeBodyEquations,init_params,time_span,args=(G,m1,m2,m3))

r1_sol=three_body_sol[:,:3]
r2_sol=three_body_sol[:,3:6]
r3_sol=three_body_sol[:,6:9]

#Erstellen der Figur
fig=plt.figure(figsize=(15,15))
#Erstellen der Achsen
ax=fig.add_subplot(111,projection="3d")
#Ploten der Orbits
ax.plot(r1_sol[:,0],r1_sol[:,1],r1_sol[:,2],color="darkblue")
ax.plot(r2_sol[:,0],r2_sol[:,1],r2_sol[:,2],color="tab:red")
ax.plot(r3_sol[:,0],r3_sol[:,1],r3_sol[:,2],color="tab:green")
#Plotten der finalen Position der Körper
ax.scatter(r1_sol[-1,0],r1_sol[-1,1],r1_sol[-1,2],color="darkblue",marker="o",s=100,label="Alpha Centauri A")
ax.scatter(r2_sol[-1,0],r2_sol[-1,1],r2_sol[-1,2],color="tab:red",marker="o",s=100,label="Alpha Centauri B")
ax.scatter(r3_sol[-1,0],r3_sol[-1,1],r3_sol[-1,2],color="tab:green",marker="o",s=100,label="Third Star")
#Hinzufügen der Beschriftungen
ax.set_xlabel("x-Koordinate",fontsize=14)
ax.set_ylabel("y-Koordinate",fontsize=14)
ax.set_zlabel("z-Kordinate",fontsize=14)
ax.set_title("Visualisierung der Orbits von Objekten im Raum\n",fontsize=14)
ax.legend(loc="upper left",fontsize=14)

ani = animation.FuncAnimation(fig, ThreeBodyEquations, frames=1000, interval=50)
plt.show()