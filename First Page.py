from tkinter import *
from PIL import ImageTk,Image
import matplotlib.pyplot as plt
import numpy as np


from matplotlib.collections import EventCollection

#####Code to Create Input Screen#####
window = Tk()

window.title("Ballistics Drag Coefficient Calculator")
window.geometry('700x500')

Header1= Label(window, text="          Range Data (Yards)")
Header1.grid(column=0, row=0)

X1lbl= Label(window, text= "X1")
X1lbl.grid(column=0, row=1)
X2lbl= Label(window, text= "X2")
X2lbl.grid(column=0, row=3)
X3lbl= Label(window, text= "X3")
X3lbl.grid(column=0, row=5)
X4lbl= Label(window, text= "X4")
X4lbl.grid(column=0, row=7)
X5lbl= Label(window, text= "X5")
X5lbl.grid(column=0, row=9)

Header2= Label(window, text="          Elevation Data (Inches)")
Header2.grid(column=3, row=0)

Y1lbl= Label(window, text= "Y1")
Y1lbl.grid(column=3, row=1)
Y2lbl= Label(window, text= "Y2")
Y2lbl.grid(column=3, row=3)
Y3lbl= Label(window, text= "Y3")
Y3lbl.grid(column=3, row=5)
Y4lbl= Label(window, text= "Y4")
Y4lbl.grid(column=3, row=7)
Y5lbl= Label(window, text= "Y5")
Y5lbl.grid(column=3, row=9)

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
txtY1.grid(column=4, row=1)
txtY2= Entry(window,width=5)
txtY2.grid(column=4, row=3)
txtY3= Entry(window,width=5)
txtY3.grid(column=4, row=5)
txtY4= Entry(window,width=5)
txtY4.grid(column=4, row=7)
txtY5= Entry(window,width=5)
txtY5.grid(column=4, row=9)

#######################################

#####Create Data entry Window For Muzzle Velocity#####

Header1= Label(window, text="          Muzzle Velocity (Ft/sec)")
Header1.grid(column=5, row=0)

Volbl= Label(window, text= "Vo")
Volbl.grid(column=5, row=1)

Muzzl= Entry(window,width=5)
Muzzl.grid(column=6, row=1)
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

BangBang = Button(window, text="Bang! Bang!", command=ClickBang)
BangBang.grid(column=0, row=11)


#####################################

#####Code To Plot Results############

x1=0

x = [x1, x2, x3, x4, x5]
y = [y1, y2, y3, y4, y5]

plt.plot(x,y, linestyle='', marker='o', color='b')
plt.xlabel('Distance in Yards')
plt.ylabel('Bullet Path in Inches')

plt.title('Ballistic Data for Brenneke K.O. Slugs')
plt.show()

####################################
window.mainloop()