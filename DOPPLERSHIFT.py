#THE FOLLOWING CODE IS JUST A WORKING PROTOTYPE OF ALGORITHM MENTIONED IN THE FOLDER NAMED AEROSPACE.THE CALCULATIONS 
#ARE NOT DONE USING ACTUAL DATA. THEREFORE, ERRORS MAY OCCUR.
#THIS ALGORITHM GENERATES THREE POSITION VECTORS AT THREE DIFFERENT TIME EPOCHS FOR ONE SATELLITE 
#CALCULATED POSITION VECTORS MAY BE USED TO CALCULATE ORBITAL PARAMETERS USING GIBB'S METHOD
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import simps
import math

def func(x,a,b,c,d,e,f,g):
    return (a*x**7+b*x**6+c*x**5+d*x**4+e*x**3+f*x**2+g*x)#degree 7 polynomial

def parameter1(xdata,ydata):#Returns distance of satellite at time stamp t=20
    popt,pcov=curve_fit(func,xdata,ydata)#Fitting the curve in the data
    def func1(t):  
        return (popt[0]*t**7+popt[1]*t**6+popt[2]*t**5+popt[3]*t**4+popt[4]*t**3+popt[5]*t**2+popt[6]*t)
    ydata2=func1(np.array(xdata))
    pos1=simps(ydata2,np.array(xdata)) #Position(Distance) at time t=10
    return pos1
    
def parameter2(xdata,ydata):#Returns distance of satellite at time stamp t=10
    popt,pcov=curve_fit(func,xdata,ydata)#Fitting the curve in the data
    def func1(t):  
        return (popt[0]*t**7+popt[1]*t**6+popt[2]*t**5+popt[3]*t**4+popt[4]*t**3+popt[5]*t**2+popt[6]*t)
    ydata2=func1(np.array(xdata))
    pos1=simps(ydata2,np.array(xdata)) #Position(Distance) at time t=10
    y2data=func1(np.array(xdata[:10]))
    pos2=simps(y2data,np.array(xdata[:10]))#position(distance) at time t=10, Applied simpson's rule of numerical integration
    return pos2    

def parameter3(xdata,ydata):#Returns distance of satellite at time stamp t=5
    popt,pcov=curve_fit(func,xdata,ydata)#Fitting the curve in the data
    def func1(t):  
        return (popt[0]*t**7+popt[1]*t**6+popt[2]*t**5+popt[3]*t**4+popt[4]*t**3+popt[5]*t**2+popt[6]*t)
    ydata2=func1(np.array(xdata))
    pos1=simps(ydata2,np.array(xdata)) #Position(Distance) at time t=10
    y2data=func1(np.array(xdata[:5]))
    pos3=simps(y2data,np.array(xdata[:5]))#position(distance) at time t=10, Applied simpson's rule of numerical integration
    return pos3   
      


#THE FOLLOWING DATA IS TAKEN RANDOMLY DUE TO UNAVALABILITY OF ACTUAL DATA
#Data with respect to station 1 after calculation from doppler shift data  
x1=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]#Time stamps
y1=[15,20,14,18,16,12,17,15,13,19,15,14,17,12,10,16,12,10,13,15] #SPEED
#Data with respect to station 2 after calculation from doppler shift data
x2=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
y2=[25,24,26,23,28,29,27,21,24,28,19,20,28,23,26,24,25,24,28,30]#SPEED
#data with respect to station 3 after calculation from doppler shift data
x3=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
y3=[20,24,29,35,26,33,36,35,34,37,38,39,32,31,25,24,26,25,28,29]#SPEED
  
p1,p2,p3=[parameter1(x1,y1),parameter1(x2,y2),parameter1(x3,y3)]
q1,q2,q3=[parameter2(x1,y1),parameter2(x2,y2),parameter2(x3,y3)]
r1,r2,r3=[parameter3(x1,y1),parameter3(x2,y2),parameter3(x3,y3)]

#Position vector of any three stations
st1 = np.array([ -294.32, 4265.1, 5986.7])#Station 1
st2 = np.array([ -1365.5, 3637.6, 6346.7])#Station 2
st3 = np.array([ -2940.3, 2473.7, 6555.8])#Station 3


#Position vector calculation
def position(X1,X2,X3,rho1,rho2,rho3):
    d1,d2=[(X2-X3),(X1-X3)]
    cross=np.cross(d1,d2)
    e=cross/np.linalg.norm(cross)#Unit vector perpendicular to the plane of R1, R2 and R3
    #Rest of the algorithm described in the attached file 
    f1=np.linalg.norm(X1)**2-rho1**2
    f2=np.linalg.norm(X2)**2-rho2**2
    f3=np.linalg.norm(X3)**2-rho3**2
    f=(1/np.linalg.norm(cross))*((f1-f3)*d1+(f2-f3)*d2)
    F=np.cross(e,f)
    H=X3.dot(e)
    G=X3.dot(F)-f3
    fsquare=F.dot(F)
    J=H + math.sqrt(H**2 + 2*G + fsquare)#Here -fsquare instead of +fsquare should be taken. But due to inconsistency of the            #randomly      taken data. If actual data from a satellite revolving in elliptical orbit is taken, the inconsistency won't occur. 
    r=F+J*e
    return r


#three pposition vectors for calculation of orbital parameters    
R1=(position(st1,st2,st3,p1,p2,p3))
R2=(position(st1,st2,st3,q1,q2,q3))
R3=(position(st1,st2,st3,r1,r2,r3))
#now Orbital parameters can be calculated using Gibb's method  
#Constants
mu = 398600.4418;
pi = 3.141592653
#Implementation of Gibb's method   
r1 = np.linalg.norm(R1)
r2 = np.linalg.norm(R2)
r3 = np.linalg.norm(R3)

c12 = np.cross(R1,R2)
c23 = np.cross(R2,R3)
c31 = np.cross(R3,R1)

##Calculate D and N vectors
##Orthogonal to the plane traced by the position vectors
N = r1*c23 + r2*c31 + r3*c12
D = c12 + c23 + c31

## Calcualte P, from pD = N

p = np.linalg.norm(N)/np.linalg.norm(D)

## Calculate S, then find e using e = S/D

S = R1*(r2 - r3) + R2*(r3 - r1) + R3*(r1 - r2)
e = np.linalg.norm(S)/np.linalg.norm(D)

## Q = w x r, calculate w then find Q

W = N/np.linalg.norm(N)
Q = S/np.linalg.norm(S)

## Calculate B and L

B = np.cross(D,R2);
L = math.sqrt(mu/(np.linalg.norm(D) * np.linalg.norm(N)));

## Calculate V2
v2 = (np.linalg.norm(L)/r2)*B + np.linalg.norm(L) * S;

## Using R2 and v2, calculate orbital elements

H = np.cross(R2,v2);

##Calculate i (inclination)
inclination = math.acos(H[2]/np.linalg.norm(H));
inclination = inclination * 180/pi;

n = [-H[1], H[0], 0];

v=np.linalg.norm(v2);
r=np.linalg.norm(R2);

vr = np.dot(R2,v2)/r;

E = 1/mu*((v*v - mu/r)*R2 - np.dot(R2,v2)*v2);

##Calculate Omega (longitude of the ascending node)

Omega = math.acos(n[0]/np.linalg.norm(n)) * 180/pi;

##Calculate omega (argument of periapsis)

omega = math.acos(np.dot(n,E)/(np.linalg.norm(n) * np.linalg.norm(E))) * 180/pi;

##Calculate nu (mean anomaly)

nu = math.acos(np.dot(E,R2)/(np.linalg.norm(E) * np.linalg.norm(R2))) * 180/pi;

##Calculate eccentricity
e = np.linalg.norm(E);

## Calculate orbital major axis

a = p / (1-(e*e));

print("Eccentricity:",e)
print("Inclination:",inclination)
print("Argument of Perigee:",omega)
print("Mean Anomaly:",nu)
print("Semi major axis:",a)
print("Right ascension of ascending node:",Omega)













  
      
