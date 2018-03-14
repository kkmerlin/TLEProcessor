import numpy as np
import math

## Define constants
mu = 398600.4418;
pi = 3.141592653

##Define 3 input vectors

R1 = np.array([ -294.32, 4265.1, 5986.7])
R2 = np.array([ -1365.5, 3637.6, 6346.7])
R3 = np.array([ -2940.3, 2473.7, 6555.8])

## Calculate vector magnitudes and cross products

r1 = np.linalg.norm(R1)
r2 = np.linalg.norm(R2)
r3 = np.linalg.norm(R3)

c12 = np.cross(R1,R2)
c23 = np.cross(R2,R3)
c31 = np.cross(R3,R1)

##Calculate D and N vectors

N = r1*c23 + r2*c31 + r3*c12
D = c12 + c23 + c31

## Check for sanity


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

result =  [ e, Omega, inclination, omega, nu, a]

print "Eccentricity:",e
print "Inclination:",inclination
print "Argument of Perigee:",omega
print "Mean Anomaly:",nu
print "Semi major axis:",a
print "Right ascension of ascending node:",Omega
