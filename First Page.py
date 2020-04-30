from tkinter import *
from PIL import ImageTk,Image
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EventCollection



#####Code to Create Input Screen#####
window = Tk()

window.title("Ballistics Drag Coefficient Calculator")
window.geometry('600x350')

Header1= Label(window, text="          Range Data (Yards)")
Header1.grid(column=0, row=0)

X1lbl= Label(window, text= "X1")
X1lbl.grid(column=1, row=1)
X2lbl= Label(window, text= "X2")
X2lbl.grid(column=1, row=3)
X3lbl= Label(window, text= "X3")
X3lbl.grid(column=1, row=5)
X4lbl= Label(window, text= "X4")
X4lbl.grid(column=1, row=7)
X5lbl= Label(window, text= "X5")
X5lbl.grid(column=1, row=9)

Header2= Label(window, text="          Elevation Data (Inches)")
Header2.grid(column=3, row=0)

Y1lbl= Label(window, text= "Y1")
Y1lbl.grid(column=4, row=1)
Y2lbl= Label(window, text= "Y2")
Y2lbl.grid(column=4, row=3)
Y3lbl= Label(window, text= "Y3")
Y3lbl.grid(column=4, row=5)
Y4lbl= Label(window, text= "Y4")
Y4lbl.grid(column=4, row=7)
Y5lbl= Label(window, text= "Y5")
Y5lbl.grid(column=4, row=9)

#####Create Entry Windows For X Coordinates#####

txtX1= Entry(window,width=5)
txtX1.grid(column=2, row=1)
txtX2= Entry(window,width=5)
txtX2.grid(column=2, row=3)
txtX3= Entry(window,width=5)
txtX3.grid(column=2, row=5)
txtX4= Entry(window,width=5)
txtX4.grid(column=2, row=7)
txtX5= Entry(window,width=5)
txtX5.grid(column=2, row=9)

###############################################

#####Create Entry Windows For Y Coordinates#####

txtY1= Entry(window,width=5)
txtY1.grid(column=5, row=1)
txtY2= Entry(window,width=5)
txtY2.grid(column=5, row=3)
txtY3= Entry(window,width=5)
txtY3.grid(column=5, row=5)
txtY4= Entry(window,width=5)
txtY4.grid(column=5, row=7)
txtY5= Entry(window,width=5)
txtY5.grid(column=5, row=9)

#######################################

#####Create Data entry Window For Muzzle Velocity#####

Header1= Label(window, text="          Muzzle Velocity (Ft/sec)")
Header1.grid(column=0, row=10)

Volbl= Label(window, text= "Vo")
Volbl.grid(column=1, row=11)

Muzzl= Entry(window,width=5)
Muzzl.grid(column=2, row=11)
######################################################

#####Create Data Entry Window For Projectile Mass#####

Header1= Label(window, text="          Projectile Mass (grains)")
Header1.grid(column=0, row=12)

masslbl= Label(window, text= "Mass")
masslbl.grid(column=1, row=13)

mass= Entry(window,width=5)
mass.grid(column=2, row=13)

######################################################

#####Create Data Entry Window For Caliber #####

Header1= Label(window, text="          Caliber (inch)")
Header1.grid(column=0, row=14)

callbl= Label(window, text= "Cal")
callbl.grid(column=1, row=15)

cal= Entry(window,width=5)
cal.grid(column=2, row=15)

######################################################

#####Create Data Entry Window For Temperature#####

Header1= Label(window, text="          Temperature (deg F)")
Header1.grid(column=3, row=10)

Templb= Label(window, text= "Temp")
Templb.grid(column=4, row=11)

Temp= Entry(window,width=5)
Temp.grid(column=5, row=11)

######################################################

#####Create Data Entry Window For Humidity#####

#Header1= Label(window, text="          Humidity (%)")
#Header1.grid(column=3, row=12)

#Humidlb= Label(window, text= "%")
#Humidlb.grid(column=4, row=13)

#Humid= Entry(window,width=5)
#Humid.grid(column=5, row=13)

######################################################

#####Create Data Entry Window For Atmospheric Pressure#####

Header1= Label(window, text="          Atmospheric Pressure (in Hg)")
Header1.grid(column=3, row=14)

Presslbl= Label(window, text= "in Hg")
Presslbl.grid(column=4, row=15)

Press= Entry(window,width=5)
Press.grid(column=5, row=15)

######################################################

#####Create Button to Enter Data and Execute Calculations#####

def ClickBang():

    X1= txtX1.get()
    x1= float(X1)

    X2= txtX2.get()
    x2= float(X2)

    X3= txtX3.get()
    x3= float(X3)

    X4= txtX4.get()
    x4= float(X4)

    X5= txtX5.get()
    x5= float(X5)

    Y1= txtY1.get()
    y1= float(Y1)

    Y2= txtY2.get()
    y2= float(Y2)

    Y3= txtY3.get()
    y3= float(Y3)

    Y4= txtY4.get()
    y4= float(Y4)

    Y5= txtY5.get()
    y5= float(Y5)

#####Code To Define Input Variables and do Unit Conversions#####
    Vel= Muzzl.get()
    Velocity= float(Vel) * .3048 #Convert to m/s

    tfin = x5 / Velocity

    Cd= float(1) #First guess of drag Coefficient

    M= mass.get()
    m= float(M)*.000064799 #Convert Projectile Mass To Kg

    D= cal.get()
    d= float(D) #Obtains diameter in inches
    A= (3.14159265359 * ((d/2)**2)) *.00064516 #Cross Sectional Area of Projectile in m^2

    g=9.807 #m/s^2

    PT= Press.get()
    pt= float(PT)* 3386.39 #converts pressure from in Hg to Pascal

    t= Temp.get()
    T= (float(t)+459.67)*(5/9) #Converts Temperature from F to K

    #####Calculates Density of Dry Air#####

    ro=pt/(286.9*T)

    #######################################

    ##ro=(PA/(287.05*T))+(PW/(461.495*T))
    ##ro=1.609


    #####Allocation of space for Derivative Estimates for Newton's Method#####
    result = np.zeros(5)
    error = np.zeros(5)
    esq = np.zeros(5)
    SESQ3= float(0)
    SESQ2= float(0)
    SESQ1= float(0)

    ##########################################################################

    #####Estimate Time Of Flight, Establish matricies for calculations#####


    h=.0002
    N= int((1.5* tfin) / h)

    #####Preallocate Space For Coordinate And Velocity Data#####

    xx = np.zeros(N)
    yy = np.zeros(N)
    v = np.zeros(N)
    a = np.zeros(N)
    vx = np.zeros(N)
    vy = np.zeros(N)
    dx = np.zeros(N)
    dy = np.zeros(N)
    dv = np.zeros(N)

    ###########################################################
    for q in range(0,25):
        #####Great Big For Loop#####
        for p in range(0, 4):

            #####Set First Values in Arrays#####
            k= 0.5*ro*Cd*A*(1/m)

            xx[0]= x1*0.9144
            yy[0]= y1*0.0254
            v[0] = Velocity

            if p == 0:
                if q == 0:
                    a[0] = ((y2-y1)/(x2-x1))*(0.0277778) #makes sure y in inches becomes yards like x

            vx[0] = np.cos(a[0])*Velocity
            vy[0] = np.sin(a[0])*Velocity

            dv[0] = (v[0]**2)*k

            dx[0] = -np.cos(a[0])*dv[0]
            dy[0] = -np.sin(a[0])*dv[0]

            ####################################

            #####Loop Through Matrix For Updated Velocity and Coordinate Data#####

            for i in range(0,N-1):

                vx[i+1]= vx[i]+(h*dx[i])
                vy[[i+1]]= vy[i]+(h*dy[i])
                v[i+1]= np.sqrt((vx[i+1]**2)+(vy[i+1]**2))
                dv[i+1]= (v[i+1]**2)*k
                a[i+1]= np.arctan(vy[i+1]/vx[i+1])
                dx[i+1]= -np.cos(a[i+1])*dv[i+1]
                dy[i+1]= (-np.sin(a[i+1])*dv[i+1])-9.807

            for n in range(0, N-1):
                xx[n + 1] = (((2/3) * (xx[n] + vx[n] * h)) + ((1/3) * (xx[n] + (h/2) * (vx[n] + vx[n+1]))));
                yy[n + 1] = (((2/3) * (yy[n] + vy[n] * h)) + ((1/3) * (yy[n] + (h/2) * (vy[n] + vy[n+1]))));

            ######################################################################

            ################################################################

            #####Code To Plot Input Data############

            x = [x1, x2, x3, x4, x5]
            y = [y1, y2, y3, y4, y5]
            xxx=xx*1.09361
            yyy=yy*39.37

            DCD=0.002

            for l in range (0,4):
                result[l] = np.interp(x[l],xxx,yyy)
                error[l] = y[l]-result[l]
                esq[l]= error[l]**2


            if p==0:
                SESQ2=sum(esq)
                Cd= Cd-DCD

            if p==1:
                SESQ1=sum(esq)
                Cd= Cd+(2*DCD)

            if p==2:
                SESQ3=sum(esq)
                FP=(SESQ3-SESQ1)/(2*DCD)
                FPP = (SESQ3 - 2 * SESQ2 + SESQ1) / (DCD**2)
                Cmod = FP / FPP
                Cd = Cd - DCD - Cmod

            if p==3:
                a[0]=a[0]+((y2-result[1])*0.0277778)/(x2-x1)


    plt.plot(x, y, linestyle='', marker='o', color='b')
    plt.plot(xxx,yyy)
    plt.xlabel('Distance in Yards')
    plt.ylabel('Bullet Path in Inches')
    CD= str(round(Cd, 4))
    TITLE = ('Drag Coefficient Cd calculated to be '+ CD)
    plt.title(TITLE)
    plt.show()


####################################

#####Code To Plot Final Curve From Model#####

#############################################

BangBang = Button(window, text="Bang! Bang!", command=ClickBang)
BangBang.grid(column=0, row=16)


#####################################


window.mainloop()
