# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 16:02:20 2019

@author: usuario
"""

import numpy as np
from scipy.optimize import fsolve
import scipy.integrate
import matplotlib.pyplot as plt


#ac=9.8 #Acceleració característica
#xc=1000 #Posició característica NO S'UTILITZA
#tc=1 #Temps característic

N=100 #Nombre d'intervals al qual es divideix cada interval de temps h

pi = np.pi

#Llegim els valors del fitxer csv                Notar que els eixos del mòbil són diferents als eixos de referència      Utilitzar a o g's?
ay = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=4) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició
az = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=5) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició
ax = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=6) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició
wy = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=7)
wz = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=8)
wx = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=9)
t = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=0) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició


#Definim el nombre  total de dades que s'han pres
line_count = len(t) #Podríem haver posat ax, ay, az o t

#Definim les arrays de velocitats, posicions i angles de desviació dels eixos. Considerem les velocitats inicials zero    
vx=np.zeros(line_count)
vy=np.zeros(line_count)
vz=np.zeros(line_count)

x=np.zeros(line_count)
y=np.zeros(line_count)
z=np.zeros(line_count)

ox=np.zeros(line_count)
oy=np.zeros(line_count)
oz=np.zeros(line_count)

a=np.array([ax, ay, az]) #Matriu on hi ha les accleracions. A cada columna hi ha les acceleracions a un instant de temps
w=np.array([wx, wy, wz]) #Matriu on hi ha les velocitats angulars. A cada columa hi ha les velocitats angulars a un instant de temps

#Passem els arrays anteriors a matrius
a = np.mat(a)
w = np.mat(w)

nax=np.zeros(line_count) #Arrays on hi haurà les acceleracions en el sistema de referència inicial
nay=np.zeros(line_count)
naz=np.zeros(line_count)

nwx=np.zeros(line_count) #Arrays on hi haurà les velocitats angulars en el sistema de referència inicial
nwy=np.zeros(line_count)
nwz=np.zeros(line_count)
    
j=0  #Definim un comptador de control


file=open("Position_data.txt", "w+") #Creem el fitxer on hi escriurem les coordenades i velocitats
file.write("#Temps  x(m)  y(m)  z(m)  vx(m/s)  vy(m/s)  vz(m/s)  ax(m/s2)  ay(m/s2)  az(m/s2)  ox(rad)  oy(rad)  oz(rad)\n")
           
           
  ##########################################################################################################
  #####################################Càlcul amb mètode propi #############################################
  ##########################################################################################################

#Definim les matrius de rotació
def rotx(thetax):
    return np.array([[1, 0, 0], [0, np.cos(thetax), -np.sin(thetax)], [0, np.sin(thetax), np.cos(thetax)]])
def roty(thetay):
    return np.array([[np.cos(thetay), 0, np.sin(thetay)], [0, 1, 0], [-np.sin(thetay), 0, np.cos(thetay)]])
def rotz(thetaz):
    return np.array([[np.cos(thetaz), -np.sin(thetaz), 0], [np.sin(thetaz), np.cos(thetaz), 0], [0, 0, 1]])

def rot(thetax, thetay, thetaz): #Definim la matriu de rotació genèrica. g=ox, b=oy, a=oz
    rotx = np.array([[1, 0, 0], [0, np.cos(thetax), -np.sin(thetax)], [0, np.sin(thetax), np.cos(thetax)]])
    roty = np.array([[np.cos(thetay), 0, np.sin(thetay)], [0, 1, 0], [-np.sin(thetay), 0, np.cos(thetay)]])
    rotz = np.array([[np.cos(thetaz), -np.sin(thetaz), 0], [np.sin(thetaz), np.cos(thetaz), 0], [0, 0, 1]])

    #Passem les arrays a matrius. Això ens permet fer productes matricials
    rotxm = np.mat(rotx)
    rotym = np.mat(roty)
    rotzm = np.mat(rotz)

    return rotzm*rotym*rotxm #Producte matricial de les matrius de rotació dels eixos. Retorna la matriu de rotació global

    #return np.array([[np.cos(a)*np.cos(b), np.cos(a)*np.sin(b)*np.sin(g)-np.sin(a)*np.cos(g), np.cos(a)*np.sin(b)*np.cos(g)+np.sin(a)*np.sin(g)], [np.sin(a)*np.cos(b), np.sin(a)*np.sin(b)*np.sin(g)+np.cos(a)*np.cos(g), np.sin(a)*np.sin(b)*np.cos(g)-np.cos(a)*np.sin(g)], [-np.sin(b), np.cos(b)*np.sin(g), np.cos(b)*np.cos(g)]])

def desrot(k):
    return np.transpose(k)


for i in range(line_count): #Trobem les desviacions dels eixos respecte el sistema de referència inicial
    if i==0:
          ox[i]=0
          oy[i]=0
          oz[i]=0
     
    else:
        #Per cada instant de temps resolem un sistema d'equacions per trobar els valors de la velocitat angular i l'angle al sistema de referència definit per 
          def equations(p):
              wx, wy, wz, anglex, angley, anglez = p
              f1 = wx-np.float(np.dot(rot(-anglex, -angley, -anglez), w[:,i])[0])
              f2 = wy-np.float(np.dot(rot(-anglex, -angley, -anglez), w[:,i])[1])
              f3 = wz-np.float(np.dot(rot(-anglex, -angley, -anglez), w[:,i])[2])
              f4 = anglex - ox[i-1] - np.float(np.dot(rot(-anglex, -angley, -anglez), w[:,i])[0]) * w[0, i]* (t[i]-t[i-1])
              f5 = angley - oy[i-1] - np.float(np.dot(rot(-anglex, -angley, -anglez), w[:,i])[1]) * w[1, i]* (t[i]-t[i-1])
              f6 = anglez - oz[i-1] - np.float(np.dot(rot(-anglex, -angley, -anglez), w[:,i])[2]) * w[2, i]* (t[i]-t[i-1])
          
              return (f1, f2, f3, f4, f5, f6)
          
          nwx[i], nwy[i], nwz[i], ox[i], oy[i], oz[i] = fsolve(equations, (0, 0, 0, 0, 0, 0))
          
         
          if ox[i]>2*pi: #Corregim els valors pels angles perquè estiguin entre -2pi i 2pi rad. Si estem parlant d'angles això no hauria d'afectar a càlculs posteriors amb aquestes variables
              ox[i] = ox[i]-2*pi
          if oy[i]>2*pi:
              oy[i] = oy[i]-2*pi
          if oz[i]>2*pi:
              oz[i] = oz[i]-2*pi
          if ox[i]<-2*pi:
              ox[i] = ox[i]+2*pi
          if oy[i]<-2*pi:
              oy[i] = oy[i]+2*pi
          if oz[i]<-2*pi:
              oz[i] = oz[i]+2*pi
          #Pels valors finals de ox, oy, oz sembla que els reusltats són bons
                
        #Rotem les acceleracions   
     
    nac = np.dot(rot(-ox[i], -oy[i], -oz[i]), a[:,i])
    
    #Escrivim els resultats de les acceleracions rotades (és a dir, totes respecte el mateix sistema de referència) en arrays per a que puguin ser accessibles al següent loop
    nax[i] = float(nac[0])
    nay[i] = float(nac[1])
    naz[i] = float(nac[2])
    
    #print(np.sqrt(a[0][i]**2+a[1][i]**2+a[2][i]**2), np.sqrt(nax[i]**2+nay[i]**2+naz[i]**2))
    #Com que a la crida incial hem redefinit els eixos, els angles són respecte el sistema de referència inical.
'''
for i in range(line_count): #Tenim els angles dins les respectives arrays ox oy oz
        if i==0: #Definim els valors incials
            x[i]=0
            y[i]=0
            z[i]=0
            
            vx[i]=0
            vy[i]=0
            vz[i]=0
      
        else:

            #Trobem la nova velocitat
            vx[i] = vx[i-1] + (t[i]-t[i-1])*(nax[i]-nax[i-1])/2
            vy[i] = vy[i-1] + (t[i]-t[i-1])*(nay[i]-nay[i-1])/2
            vz[i] = vz[i-1] + (t[i]-t[i-1])*(naz[i]-naz[i-1])/2
            #vx[i] = vx[i-1] + (t[i]-t[i-1])*nax[i]
            #vy[i] = vy[i-1] + (t[i]-t[i-1])*nay[i]
            #vz[i] = vz[i-1] + (t[i]-t[i-1])*naz[i]
            
#            vx[i] = vx[i-1] + (t[i]-t[i-1])*(nax[i]-nax[i-1])/2 + nax[i-1]*(t[i]-t[i-1])
#            vy[i] = vy[i-1] + (t[i]-t[i-1])*(nay[i]-nay[i-1])/2 + nay[i-1]*(t[i]-t[i-1])
#            vz[i] = vz[i-1] + (t[i]-t[i-1])*(naz[i]-naz[i-1])/2 + naz[i-1]*(t[i]-t[i-1])
            
#            vx[i] = vx[i-1] + (t[i]-t[i-1])*(ax[i]-ax[i-1])/2 + ax[i-1]*(t[i]-t[i-1])
#            vy[i] = vy[i-1] + (t[i]-t[i-1])*(ay[i]-ay[i-1])/2 + ay[i-1]*(t[i]-t[i-1])
#            vz[i] = vz[i-1] + (t[i]-t[i-1])*(az[i]-az[i-1])/2 + az[i-1]*(t[i]-t[i-1])

        #Primer rotar, després moure's
            #Situem les posicions en un vector de velocitat i hi apliquem la rotació
 #           v = np.array([nvx, nvy, nvz])  #CAL FER EL CÀLCUL HAVENT ROTAT LA VELOCITAT, AMB AIXÒ JA ENS SORTIRÀ LA POSICIÓ ROTADA  TAMBÉ
  #          Rot = np.dot(rotx(ox[i]), np.dot(roty(oy[i]), rotz(oz[i]))) #Matriu de rotació
   #         nvel = np.dot(Rot, v)  #Vector posició amb la rotació

    #        vx[i] = nvel[0]
     #       vy[i] = nvel[1]
      #      vz[i] = nvel[2]
            
#           #Amb la nova velocitat, trobem la nova posició
            
#            x[i] = x[i-1] + (t[i]-t[i-1])*(vx[i]-vx[i-1])/2
#            y[i] = y[i-1] + (t[i]-t[i-1])*(vy[i]-vy[i-1])/2 
#            z[i] = z[i-1] + (t[i]-t[i-1])*(vz[i]-vz[i-1])/2

#            x[i] = (t[i]-t[i-1])*(vx[i]-vx[i-1])/2 + vx[i-1]*(t[i]-t[i-1])
#            y[i] = (t[i]-t[i-1])*(vy[i]-vy[i-1])/2 + vy[i-1]*(t[i]-t[i-1])
#            z[i] = (t[i]-t[i-1])*(vz[i]-vz[i-1])/2 + vz[i-1]*(t[i]-t[i-1])

            #x[i] = x[i-1] + (t[i]-t[i-1])*(vx[i]-vx[i-1])/2 + vx[i-1]*(t[i]-t[i-1])
            #y[i] = y[i-1] + (t[i]-t[i-1])*(vy[i]-vy[i-1])/2 + vy[i-1]*(t[i]-t[i-1])
            #z[i] = z[i-1] + (t[i]-t[i-1])*(vz[i]-vz[i-1])/2 + vz[i-1]*(t[i]-t[i-1])
          
            #Assignem els nous valors a les posicions que es posen d'output
           # x[i] = npos[0]
           # y[i] = npos[1]
           # z[i] = npos[2]
            
            #if j%100==0:
             #   print(x[i], y[i], z[i], j)
#                print("mod(vx) %f , %d\n" %(np.sqrt(vx[i]**2+vy[i]**2+vz[i]**2), j))
                       
            #vx[i] = vx[i-1] + nax[i]*(t[i]-t[i-1])
            #vy[i] = vy[i-1] + nay[i]*(t[i]-t[i-1])
            #vz[i] = vz[i-1] + naz[i]*(t[i]-t[i-1])
            x[i] = x[i-1] + vx[i-1]*(t[i]-t[i-1]) + nax[i] * (t[i]-t[i-1])**2
            y[i] = y[i-1] + vy[i-1]*(t[i]-t[i-1]) + nay[i] * (t[i]-t[i-1])**2
            z[i] = z[i-1] + vz[i-1]*(t[i]-t[i-1]) + naz[i] * (t[i]-t[i-1])**2
                         
        j += 1               
        #file.write("%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f \n" %(t[i],x[i],y[i],z[i],vx[i], vy[i], vz[i], ax[i],ay[i],az[i], ox[i], oy[i], oz[i]))
'''
vx = scipy.integrate.cumtrapz(nax, x=t, initial=0.)
vy = scipy.integrate.cumtrapz(nay, x=t, initial=0.)
vz = scipy.integrate.cumtrapz(naz, x=t, initial=0.)
x = scipy.integrate.cumtrapz(vx, x=t, initial=0.)
y = scipy.integrate.cumtrapz(vy, x=t, initial=0.)
z = scipy.integrate.cumtrapz(vz, x=t, initial=0.)

#print(ox[line_count-1], oy[line_count-1], oz[line_count-1]) 
print(a)
print(nac)
 
file.close()
     
print("El fitxer té ", j, " punts.")

fig, axis = plt.subplots(4, 1, figsize=(15,12), sharex='all')

axis[0].plot(t, nax)
axis[0].plot(t, nay)
axis[0].plot(t, naz)
axis[1].plot(t, vx*3.6)
axis[1].plot(t, vy*3.6)
axis[1].plot(t, vz*3.6)
axis[2].plot(t, wx)
axis[2].plot(t, wy)
axis[2].plot(t, wz)
axis[3].plot(t, ox)
axis[3].plot(t, oy)
axis[3].plot(t, oz)

axis[0].set_ylabel("Lin. acc. (m/s2)")
axis[1].set_ylabel("Lin. vel. (km/h)")
axis[2].set_ylabel("Ang. vel. (rad/s)")
axis[3].set_ylabel("Angle (rad)")
axis[3].set_xlabel("Time (s)")

axis[3].set_yticks([-2*pi,-pi, 0, pi, 2*pi])
axis[3].set_yticklabels([str(r'-2$\pi$'), str(r'$-\pi$'),str(r'0'), str(r'$\pi$') ,str(r'2$\pi$')])


for i in range(4):
    axis[i].grid()
    
plt.savefig("Profiles.jpeg")

'''
ytick_angle = np.arange(-2*np.pi, 2*np.pi+0.5*np.pi, 0.5*np.pi)
ylabel_angle = [r"$" + format(r, ".2g")+ r"\pi$" for r in ytick_angle]
axis[3].set_yticks(ytick_angle*np.pi)
axis[3].set_yticklabels(ylabel_angle, fontsize=20)

fig_pos = plt.figure(2)
ax_pos = plt.axes(projection='3d')
ax_pos.plot3D(x, y, z)
'''

plt.show()

            
#            ##############################################################################################
#    ###################################### Càlcul amb mètode RK4 #########################
#    ###############################################################################################
#
#          
##Càlcul amb mètode RK4
#    
##Normalitzem els valors -
#
#ax = ax/ac
#ay = ay/ac
#az = az/ac
#t = t/tc
#
##Calculem
#for i in range(line_count):
#              
#    if i==0: #Escrivim les coordenades inicials
#        nx=0 
#        ny=0
#        nz=0
#        nvx=0
#        nvy=0
#        nvz=0
#        
#    else:
#        #Definim l'interval entre temps h
#        h=t[i]-t[i-1]
#    
#        #Definim els coeficients       
#        
#        c0x=vx[i]           #No sé si està bé perquè vx[i]=0 en tot moment. No estic agafant el punt anterior
#        k0x=ax[i-1]             #PROVAR MÈTODE HEUN
#        c0y=vy[i]
#        k0y=ay[i-1]
#        c0z=vz[i]
#        k0z=az[i-1]
#            
#        c1x=vx[i]+h*k0x/2
#        k1x=ax[i-1]
#        c1y=vy[i]+k0y*h/2
#        k1y=ay[i-1]
#        c1z=vz[i]+k0z*h/2
#        k1z=az[i-1]
#        
#        c2x=vx[i]+k1x*h/2
#        k2x=ax[i-1]
#        c2y=vy[i]+k1y*h/2
#        k2y=ay[i-1]
#        c2z=vz[i]+k1z*h/2
#        k2z=az[i-1]
#        
#        c3x=vx[i]+h*k2x
#        k3x=ax[i-1]
#        c3y=vy[i]+h*k2y
#        k3y=ay[i-1]
#        c3z=vz[i]+h*k2z
#        k3z=az[i-1]
#        
#        nx=x[i-1]+(c0x + 2*c1x + 2*c2x + c3x)*h/6
#        nvx=vx[i-1]+(k0x + 2*k1x + 2*k2x + k3x)*h/6
#        ny=y[i-1]+(c0y + 2*c1y + 2*c2y + c3y)*h/6
#        nvy=vy[i-1]+(k0y + 2*k1y + 2*k2y + k3y)*h/6
#        nz=z[i-1]+(c0z  +2*c1z + 2*c2z + c3z)*h/6
#        nvz=vz[i-1]+(k0z + 2*k1z + 2*k2z + k3z)*h/6
#        
#        
#    x[i]=nx
#    y[i]=ny
#    z[i]=nz
#        
#    vx[i]=nvx
#    vy[i]=nvy
#    vz[i]=nvz
#        
#    j += 1
#               
#    file.write("%f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n" %(t[i],x[i],y[i],z[i],vx[i], vy[i], vz[i], ax[i],ay[i],az[i]))
#    
#file.close()
#     
#print("El fitxer té ", j, " punts.")
#



#              ##############################################################################################
#    ###################################### Càlcul amb mètode d'Euler #########################
#    ###############################################################################################
#
#for i in range(line_count):
#    
#    if i==0:
#        x[i]=0
#        
#    elif i==1:
#        
#        h=t[i]-t[i-1]
#        x[i]=ax[i]*h**2+2*x[i-1]
#        
#    
#    else: 
#        h=t[i]-t[i-1]
#        print(h)
#        x[i]=ax[i]*h**2+2*x[i-1]-x[i-2]
#        
#                
#    i += 1
#print(line_count)
#
##Escrivim els resultats en un fitxer de dades
#for i in range(line_count):
#    file.write("%f  %f  %f\n" %(t[i],x[i], ax[i]))





    ##############################################################################################
    ###################################### Càlcul amb mètode de Bolzano #########################
    ###############################################################################################
                #NO S'HA POGUT RESOLDRE PER BOLZANO PERQUÈ NO SEMPRE LA FUNCIÓ TALLA L'EIX x

#N=450 #Nombre màxim d'iteracions per trobar les posicions
#
##Llegim els valors del fitxer csv
#ax = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=4) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició
#ay = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=5) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició
#az = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=6) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició
#t = np.loadtxt("BTM.csv", dtype=float, delimiter=',', usecols=0) #és una array on hi ha guardats els valors de la funció ona misteriosa a cada posició
#
##Definim el nombre  total de dades que s'han pres
#line_count = len(ax) #Podríem haver posat ax, ay, az o t
#
##Definim les arrays de velocitats i posicions. Considerem les velocitats inicials zero    
#vx=np.zeros(line_count)
#vy=np.zeros(line_count)
#vz=np.zeros(line_count)
#
##x=np.zeros(line_count)
##y=np.zeros(line_count)
##z=np.zeros(line_count)
#
#x=np.zeros(line_count)
#y=np.zeros(line_count)
#z=np.zeros(line_count)
#    
#j=0  #Definim un comptador de control
#
#
#file=open("Position_data.txt", "w+") #Creem el fitxer on hi escriurem les coordenades i velocitats
#file.write("#Temps  x  y  z  vx  vy  vz  ax  ay  az\n")
#
##Definim la posició que ens dona la coordenada nova        
#def coord(r, u, k): #p1=posició nova p0=posició vella c=p0**2-a*h
#    return r**2 + 2*r*u + k 
#
#d = 0.3 #Increment de posició
#e = 0.01 #error
#
##Definim un valor amb el qual compararem el signe de la funció a la primera iteració
#if coord(d, 0, d**2-ax[0]*(t[1]-t[0]))>0:
#    anscoord = 1
#
#if coord(d, 0, d**2-ax[0]*(t[1]-t[0]))<0:
#    anscoord = -1
#
#
##Càlcul amb el mètode de Bolzano
#for i in range(line_count):
#    
#    if i==0:
#        
#        x[i]=0
#        y[i]=0
#        z[i]=0
#        
#    else:    
#
#        h=t[i]-t[i-1]
#        c=x[i-1]**2-ax[i]*h**2
#        
#        print(h)
#        print(c)
#    
#        
#        nx = x[i-1]
#        print(nx)
#        
#                     #FEM EL CÀLCUL PER A X I DESPRÉS JA ESCRIUREM PER A Y,Z
#    
#        while q>e: #Això està bé?
#            
#            if anscoord*coord(nx, x[i-1], c)<0: #Canvi de signe NO S'ARRIBA A CANVIAR MAI DE SIGNE!
#                x[i]=nx 
#                break
#    
#            nx += d    
#            q=abs(coord(nx,x[i-1],c))
#            
#            
#        if coord(nx, x[i-1], d**2-ax[0]*(t[1]-t[0]))>0: #Revisar
#            anscoord = 1
#
#        if coord(nx, x[i-1], d**2-ax[0]*(t[1]-t[0]))<0:
#            anscoord = -1
#
#                
#    j += 1
#    i += 1
#              
##Escrivim els resultats en un fitxer de dades
#for i in range(line_count):
#    file.write("%f  %f  %f\n" %(t[i],x[i], ax[i]))
#    
#file.close()
#     
#print("El fitxer té ", j, " punts.")
    
   
    ##################################################################################
    ####################################################################################
    ######################################################################################

#ax=[]
#ay=[]
#az=[]
#t=[]
#x=[]
#y=[]
#z=[]

#with open('BTM.csv') as csv_file:
#    csv_reader = csv.reader(csv_file, delimiter=',')
#    line_count = 0
#    j=0    
#    
#    file=open("Position_data.txt", "w+") #Creem el fitxer on hi escriurem les coordenades
#    file.write("#Temps  x  y  z  vx  vy  vz  ax  ay  az\n")
#    for row in csv_reader:
#        line_count += 1
#        ax.append(row[4]) #Afegim els valors a la array corresponent
#        ay.append(row[5])
#        az.append(row[6])
#        t.append(row[0])
#        
#        
#    #Passem els elements llegits de str a flaot per poder operar. Canvi de índexs!
#    for i in range(1, line_count):
#        ax[i-1]=float(az[i]) #Notar que hem passat de ax[1] a ax[0]
#        ay[i-1]=float(ax[i])
#        az[i-1]=-float(ay[i])#EL PROBLEMA ÉS QUE AQUESTS VALORS NO ELS AGAFA BÉ
#        t[i-1]=float(t[i])
#        
#        ax[i]=float(az[i]) #Notar que hem passat de ax[1] a ax[0]
#        ay[i]=float(ax[i])
#        az[i]=-float(ay[i])
#        t[i]=float(t[i])
#        #Valor màxim dels índexs: line_count-1
#   
#    #Definim el temps i l'acceleració característics
#    print(ax[line_count-2]*2)
#    print(az[10])
##    ac = max(ax+ay+az)
##    tc = max(t) #De fet, és t[line_count-1]
#
#    #Normalitzem
##    for i in range(line_count-1):
##        ax[i] = float(ax[i]/ac)
##        ay[i] = ay[i]/ac
##        az[i] = az[i]/ac
##        t[i] = t[i]/tc
#        
#    #Definim les arrays de velocitats. Considerem les velocitats inicials zero    
#    vx=np.zeros(line_count)
#    vy=np.zeros(line_count)
#    vz=np.zeros(line_count)
#    
#    #Càlcul amb mètode RK4
#    for i in range(line_count):
#              
#        if i==0: #Escrivim les coordenades inicials
#             nx=0 
#             ny=0
#             nz=0
#             nvx=0
#             nvy=0
#             nvz=0
#            
#        else:
#            #Definim l'interval entre temps h
#            h=t[i]-t[i-1]
#        
#            #Definim els coeficients
#            c0x=vx[i]
#            k0x=ax[i] 
#            c0y=vy[i]
#            k0y=ay[i]
#            c0z=vz[i]
#            k0z=az[i]
#            
#            c1x=vx[i]+h*k0x/2
#            k1x=ax[i]
#            c1y=vy[i]+k0y*h/2
#            k1y=ay[i]
#            c1z=vz[i]+k0z*h/2
#            k1z=az[i]
#            
#            c2x=vx[i]+k1x*h/2
#            k2x=ax[i]
#            c2y=vy[i]+k1y*h/2
#            k2y=ay[i]
#            c2z=vz[i]+k1z*h/2
#            k2z=az[i]
#        
#            c3x=vx[i]+h*k2x
#            k3x=ax[i]
#            c3y=vy[i]+h*k2y
#            k3y=ay[i]
#            c3z=vz[i]+h*k2z
#            k3z=az[i]
#        
#            nx=x[i-1]+(c0x + 2*c1x + 2*c2x + c3x)*h/6
#            nvx=vx[i-1]+(k0x + 2*k1x + 2*k2x + k3x)*h/6
#            ny=y[i-1]+(c0y + 2*c1y + 2*c2y + c3y)*h/6
#            nvy=vy[i-1]+(k0y + 2*k1y + 2*k2y + k3y)*h/6
#            nz=z[i-1]+(c0z  +2*c1z + 2*c2z + c3z)*h/6
#            nvz=vz[i-1]+(k0z + 2*k1z + 2*k2z + k3z)*h/6
#            
#    
#        x.append(nx)
#        y.append(ny)
#        z.append(nz)
#        
#        vx[i]=nvx
#        vy[i]=nvy
#        vz[i]=nvz
#        
#        j += 1
#               
#        file.write("%f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n" %(t[i],x[i],y[i],z[i],vx[i], vy[i], vz[i], ax[i],ay[i],az[i]))
#    
#    file.close()
#     
#    print("El fitxer té ", j, " punts.")

    #print("x[", i, "] = ", x[i])
    
    #print(t[26522]*t[26523])
    #print(len(ax))
    #print(f"Hi ha {line_count} entrades")  
    
    
    
     #Escrivim les operacions que calculen les posicions i anem guardant els reusltats a l'array x,y,z corresponent   
#    for i in range(line_count-1): (o -2?)
#
#        if i==0:
#            nx=0 #Expressem què passa en els casos i=0, i=1
#            ny=0
#            nz=0
#            j += 1
#            
#        elif i==1:
#            h=t[i]-t[i-1]
#            nx=az[i-1]*h**2+2*x[i-1] #Expressem què passa en el cas i=1
#            ny=ax[i-1]*h**2+2*y[i-1]
#            nz=-ay[i-1]*h**2+2*z[i-1]
#            j += 1
#            
##        elif i==line_count-1:
##            nx=x[0]
##            ny=y[0]
##            nz=z[0]
#        
#        else:
#            h=t[i]-t[i-1]  #POT SER QUE EL PROBLEMA SIGUI QUE h NO ÉS CTT?
#            nx=az[i-1]*h**2+2*x[i-1]-x[i-2] #Expressem què passa en el cas genèric
#            ny=ax[i-1]*h**2+2*y[i-1]-y[i-2]
#            nz=-ay[i-1]*h**2+2*z[i-1]-z[i-2]
#            j += 1
     
     ##############################################################################
     #################################################################################
     #################################################################################
     
      #Càcul amb el mètode d'Adams-Bashford
#    for i in range(line_count):
#        if i==0: #Especifiquem que les posicions i velocitats inicials són 0
#            nx=0
#            ny=0
#            nz=0
#            nvx=0
#            nvy=0
#            nvz=0
#        
#        elif i==1:
#            
#            nvx=vx[i-1]+0.5*(t[i-1]+t[i])*ax[i-1]
#            nvy=vy[i-1]+0.5*(t[i-1]+t[i])*ay[i-1]
#            nvz=vz[i-1]+0.5*(t[i-1]+t[i])*az[i-1]
#            nx=x[i-1]+0.5*(t[i-1]+t[i])*vx[i-1]
#            ny=y[i-1]+0.5*(t[i-1]+t[i])*vy[i-1]
#            nz=z[i-1]+0.5*(t[i-1]+t[i])*vz[i-1]
#            
#        else:
#            nvx=vx[i-1]+0.5*(t[i-1]+t[i]-2*t[i-2])*ax[i-1]+0.5*(2*t[i-1]-t[i]-t[i-2])*ax[i-2]
#            nvy=vy[i-1]+0.5*(t[i-1]+t[i]-2*t[i-2])*ay[i-1]+0.5*(2*t[i-1]-t[i]-t[i-2])*ay[i-2]
#            nvz=vz[i-1]+0.5*(t[i-1]+t[i]-2*t[i-2])*az[i-1]+0.5*(2*t[i-1]-t[i]-t[i-2])*az[i-2]
#            
#            nx=x[i-1]+0.5*(t[i-1]+t[i]-2*t[i-2])*vx[i-1]+0.5*(2*t[i-1]-t[i]-t[i-2])*vx[i-2]
#            ny=y[i-1]+0.5*(t[i-1]+t[i]-2*t[i-2])*vy[i-1]+0.5*(2*t[i-1]-t[i]-t[i-2])*vy[i-2]
#            nz=z[i-1]+0.5*(t[i-1]+t[i]-2*t[i-2])*vz[i-1]+0.5*(2*t[i-1]-t[i]-t[i-2])*vz[i-2]    
      
      
      
      
#            ##############################################################################################
#    ###################################### Càlcul amb mètode propi (adaptació de trapezis) #########################
#    ###############################################################################################       
#           
##Càlcul
#
#for i in range(line_count):
#    if i==0: #Escrivim les coordenades incials
#        x[i] = 0
#        y[i] = 0
#        z[i] = 0
#        
#        vx[i] = 0
#        vy[i] = 0
#        vz[i] = 0
#        
#    #ESPECIFICAR CAS a[line_count]  
#    
#    elif i==line_count-1:
#        x[i]=x[i-1]
#        y[i]=y[i-1]
#        z[i]=z[i-1]
#        
#        vx[i]=vx[i-1]
#        vy[i]=vy[i-1]
#        vz[i]=vz[i-1]
#        
#    else:
#        #Trobem la recta y=m*x+c que uneix a[i] i a[i+1]
#        x1 = (t[i], ax[i])#Definim els dos punts els quals en volem trobar la recta que els uneix
#        x2 = (t[i+1], ax[i+1])
#        
#        coordx = np.polyfit(x1, x2, 1) #Especifiquem els dos punts i de quin grau és el polinomi al qual volem fer el fit
#                                    #Obtenim una array on coord[0]=m, coord[1]=c
#        
#        rx=np.linspace(ax[i], ax[i+1], N, endpoint=True)
#        T=np.linspace(t[i], t[i+1], N, endpoint=True)      #Adaptar N en funció de l'alçada del pic? CREC QUE CAL FER-HO I ÉS DE JUSTÍCIA
#        
#               
#        for k in range(N):
#            px[k] = coordx[0]*rx[k] + coordx[1] #Trobem els valors a la recta
#           Ax[k]=px[k]*T[k] #Trobem el valor de l'àrea sota la recta en el petit interval. És la seva velocitat en aquest interval
#  