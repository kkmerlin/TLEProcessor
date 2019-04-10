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
x1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]#Time stamps
y1=[6831.364579043655, 6740.850822038237, 6620.243605101301, 6458.014493658278, 6237.037598559635, 5931.685507415964, 5503.903984015603, 4899.39563649145, 4048.9044134429896, 2887.693965597149, 1409.8162620138917, 263.4108033539343, 1897.751744773451, 3283.5217264849252, 4343.574761925622, 5110.007716248858, 5652.925115176098, 6037.812014780742, 6313.659282761534, 6514.206671309588] #SPEED
#Data with respect to station 2 after calculation from doppler shift data
x2=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
y2=[7096.59754103314, 7037.794946110064, 6952.833025271406, 6828.0412760003355, 6639.7893058758145, 6345.4505647112865, 5865.562659033127, 5052.9027833895325, 3669.098407796368, 1512.062235767461, 1118.2141634749757, 3387.7121166054285, 4883.555054240137, 5766.906774004402, 6286.4391928767045, 6603.0417954405575, 6804.302893539096, 6937.076566288755, 7027.178082531641, 7089.441901351016]#SPEED
#data with respect to station 3 after calculation from doppler shift data
x3=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
y3=[6989.497374791167, 6893.071521201201, 6752.671857144409, 6542.838112750756, 6218.592359417079, 5698.778673394734, 4840.810945027412, 3434.1562543851946, 1339.449830451687, 1136.0813694594797, 3283.7254356193903, 4746.524662821048, 5642.129367082341, 6184.013614246464, 6521.036861560033, 6738.493986017758, 6883.641085744482, 6983.158194181431, 7052.646238007066, 7101.582873933399]#SPEED
  
p1,p2,p3=[parameter1(x1,y1),parameter1(x2,y2),parameter1(x3,y3)]
q1,q2,q3=[parameter2(x1,y1),parameter2(x2,y2),parameter2(x3,y3)]
r1,r2,r3=[parameter3(x1,y1),parameter3(x2,y2),parameter3(x3,y3)]

#Position vector of any three stations
st1 = np.array([ 454177.3804063904, 4321209.12405049, 4659454.433015735])#Station 1
st2 = np.array([ 558000.6604074897, 3970381.0041338456, 4951196.920442361])#Station 2
st3 = np.array([ 665795.5068229569, 3775913.9534281944, 5088106.834511303])#Station 3


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
    k=H**2 + 3*G - fsquare
    print(k)
    J=H + math.sqrt(abs(H**2 + 3*G +fsquare))#Here -fsquare instead of +fsquare should be taken. But due to inconsistency of the            #randomly      taken data. If actual data from a satellite revolving in elliptical orbit is taken, the inconsistency won't occur. 
    r=F+J*e
    return r


#three pposition vectors for calculation of orbital parameters    
R1=(position(st1,st2,st3,p1,p2,p3))
R2=(position(st1,st2,st3,q1,q2,q3))
R3=(position(st1,st2,st3,r1,r2,r3))
#now Orbital parameters can be calculated using Gibb's method  
#Constants:
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












  
      
