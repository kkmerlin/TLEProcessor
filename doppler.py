import math
m_earth=5.97219e24
const_G = 6.667384e-11
const_c = 2.99792458e8
radius_earth = 6371e3
radius_sat = 0
period_T = 0
rest_freq = 0
d_freq = 0
velocity_sat = 0
velocity_receiver=0
d_velocity = 0
print("GPS L1 signal = 1575.42MHz....ISS = 145.8MHz")
rest_freq=float(input("Enter frequency in MHz:"))
#Period of satellite orbits
print("Geostationary orbit=1440 minutes\nGPS orbit=720 minutes\nISS orbit=92 minutes")
period_T=float(input("Please enter the period of satellite orbit in minutes:"))
#Mathematical Calculations
period_T=period_T*60.0
radius_sat = (const_G*m_earth*math.pow(period_T,2)/(math.pow(2*math.pi,2)))**0.33333
sat_altitude = radius_sat - radius_earth
print("\nSatellite avg. altitude: " + str(round(sat_altitude)/1e3) + " km")
velocity_sat = 2*math.pi*radius_sat/period_T
print("\nSatellite velocity: " + str(round(velocity_sat)/1e3) + " km/s")
velocity_receiver = 2*math.pi*radius_earth/86400.0
print("\nReceiver station velocity: " + str(round(velocity_receiver)/1e3) + " km/s")
d_velocity = 2*math.pi*(radius_sat/period_T + radius_earth/86400.0) #Calculating delta velocity
d_freq = 1e3*rest_freq*d_velocity/const_c #Calculating frequency shift
print("\nFrequency Shift: (+/-)" + str(d_freq) + " kHz")


