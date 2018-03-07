# -*- coding: utf-8 -*-
import urllib
import pytz
import numpy as np
from datetime import timedelta,datetime
pdt = pytz.timezone('Asia/Kolkata')

# GM is standard gravitational parameter
GM = 398600.4418

def splitElement(tle):
    "Splits a two line element set into title and it's two lines with stripped lines"
    return map(lambda x: x.strip(), tle.split('\n'))

def CheckSum(line):
    """Apply the basic checksum algorithm.For each line the checksums are calculated by adding all the digits 
        of the line including the line number"""
    return sum(map(int, filter(lambda c: c >= '0' and c <= '9', line[:-1].replace('-','1')))) % 10


def checkValid(tle):
    "Check if the given format of TLE is valid or not with the help of CheckSum function"
    title, line1, line2 =  splitElement(tle)

    return line1[0] == '1' and line2[0] == '2' and \
           line1[2:7] == line2[2:7] and \
           int(line1[-1]) == CheckSum(line1) and int(line2[-1]) == CheckSum(line2)

def SciNotation(sn):
    "Specific format is 5 digits, a + or -, and 1 digit, ex: 01234-5 which is 0.01234e-5"
    return 0.00001*float(sn[5]) * 10**int(sn[6:])

def eccentricAnomalyFromMean(mean_anomaly, eccentricity, initValue,
                             maxIter=500, maxAccuracy = 0.0001):
    
    """All inputs as well as outputs are in radians"""
    mean_anomaly = mean_anomaly
    e0 = initValue
    for x in xrange(maxIter):
        e1 = e0 - (e0 - eccentricity * np.sin(e0) - mean_anomaly) / (1.0 - eccentricity * np.cos(e0))
        if (abs(e1-e0) < maxAccuracy):
            break
    return e1

def print_this(tle, Info = True, labels = True):
    "Returns commented information on a two line element"
    title, line1, line2 =  splitElement(tle)
    if not checkValid(tle):
        print "Invalid element."
        return

    satellite_number                                        = int(line1[2:7])
    class_of_satellite                                      = line1[7:8]
    designator_year                                         = int(line1[9:11])
    designator_launch_number                                = int(line1[11:14])
    designator_piece_of_launch                              = line1[14:17]
    epoch_year                                              = int(line1[18:20])
    epoch                                                   = float(line1[20:32])
    first_time_derivative_of_the_mean_motion_divided_by_two = float(line1[33:43])
    second_time_derivative_of_mean_motion_divided_by_six    = SciNotation(line1[44:52])
    bstar_drag_term                                         = SciNotation(line1[53:61])
    the_number_0                                            = float(line1[62:63])
    element_number                                          = float(line1[64:68])
    checksum1                                               = float(line1[68:69])

    satellite        = int(line2[2:7])
    inclination      = float(line2[8:16])
    right_ascension  = float(line2[17:25])
    eccentricity     = float(line2[26:33]) * 0.0000001
    argument_perigee = float(line2[34:42])
    mean_anomaly     = float(line2[43:51])
    mean_motion      = float(line2[52:63])
    revolution       = float(line2[63:68])
    checksum2        = float(line2[68:69])

    # Inferred Epoch date
    year = 2000 + epoch_year if epoch_year < 70 else 1900 + epoch_year
    epoch_date = datetime(year=year, month=1, day=1, tzinfo=pytz.utc) + timedelta(days=epoch-1) # Have to subtract one day to get correct midnight

    # Time difference of now from epoch, offset in radians
    diff = datetime.now().replace(tzinfo=pytz.utc) + timedelta(hours=8) - epoch_date # Offset for PDT
    diff_seconds = 24*60*60*diff.days + diff.seconds + 1e-6*diff.microseconds # sec
    print "Time offset: %s" % diff
    motion_per_sec = mean_motion * 2*np.pi / (24*60*60) # rad/sec
    print "Radians per second: %g" % motion_per_sec
    offset = diff_seconds * motion_per_sec #rad
    print "Offset to apply: %g" % offset
    mean_anomaly += offset * 180/np.pi % 360

    # Inferred period
    day_seconds = 24*60*60
    period = day_seconds * 1./mean_motion

    # Inferred semi-major axis (in km)
    semi_major_axis = ((period/(2*np.pi))**2 * GM)**(1./3)

    # Inferred true anomaly
    eccentric_anomaly = eccentricAnomalyFromMean(mean_anomaly * np.pi/180, eccentricity, mean_anomaly * np.pi/180)
    true_anomaly = 2*np.arctan2(np.sqrt(1+eccentricity) * np.sin(eccentric_anomaly/2.0), np.sqrt(1-eccentricity) * np.cos(eccentric_anomaly/2.0))
    # Convert to degrees
    eccentric_anomaly *= 180/np.pi
    true_anomaly *= 180/np.pi
    sat_angle = true_anomaly * np.pi/180
    satr = (semi_major_axis * (1-eccentricity**2)) / (1 + eccentricity*np.cos(sat_angle))
    satx = satr * np.cos(sat_angle)
    saty = satr * np.sin(sat_angle)
    satz = 0
    #calculation of rotation matrix for inclination of the satellite
    inc = inclination * np.pi / 180.;
    R = np.matrix([[1, 0, 0],
                   [0, np.cos(inc), -np.sin(inc)],
                   [0, np.sin(inc), np.cos(inc)]    ])
    #Calculation of rotation matrix for summation of argument of perigee and right ascension
    rot = (right_ascension + argument_perigee) * np.pi/180
    R2 = np.matrix([[np.cos(rot), -np.sin(rot), 0],
                    [np.sin(rot), np.cos(rot), 0],
                    [0, 0, 1]    ]) 
    #Calculation of latitude and longitude   
    sat = (R * R2 * np.matrix([satx, saty, satz]).T ).flatten()
    satx = sat[0,0]
    saty = sat[0,1]
    satz = sat[0,2]

    c = np.sqrt(satx*satx + saty*saty)
    lat = np.arctan2(satz, c) * 180/np.pi
    lon = np.arctan2(saty, satx) * 180/np.pi

    if (Info):
        print "----------------------------------------------------------------------------------------"
        print tle
        print "---"
        print "Satellite_Name                                            = %s" % title
        print "Satellite_number                                          = %g (%s)" % (satellite_number, "Unclassified" if class_of_satellite == 'U' else "Classified")
        print "International_Designator                                  = YEAR: %02d, LAUNCH #%d, PIECE: %s" % (designator_year,designator_launch_number,designator_piece_of_launch)
        print "Epoch_Date                                                = %s  (YEAR:%02d DAY:%.11g)" % (epoch_date.strftime("%Y-%m-%d %H:%M:%S.%f %Z"), epoch_year, epoch)
        print "First Time Derivative of the Mean Motion divided by two   = %g" % first_time_derivative_of_the_mean_motion_divided_by_two
        print "Second Time Derivative of Mean Motion divided by six      = %g" % second_time_derivative_of_mean_motion_divided_by_six
        print "BSTAR drag term                                           = %g" % bstar_drag_term
        print "The number 0                                              = %g" % the_number_0
        print "Element number                                            = %g" % element_number
        print
        print "Inclination [in Degrees]                                     = %g°" % inclination
        print "Right Ascension of the Ascending Node [in Degrees]           = %g°" % right_ascension
        print "Eccentricity                                              = %g" % eccentricity
        print "Argument of Perigee [in Degrees]                             = %g°" % argument_perigee
        print "Mean Anomaly [in Degrees] Anomaly                            = %g°" % mean_anomaly
        print "Eccentric_Anomaly                                         = %g°" % eccentric_anomaly
        print "True_Anomaly                                              = %g°" % true_anomaly
        print "Mean Motion [Revs per day] Motion                         = %g" % mean_motion
        print "Period                                                    = %s" % timedelta(seconds=period)
        print "Revolution number at epoch [Revs]                         = %g" % revolution

        print
        print "semi_major_axis = %gkm" % semi_major_axis
        print "eccentricity    = %g" % eccentricity
        print "inclination     = %g°" % inclination
        print "arg_perigee     = %g°" % argument_perigee
        print "right_ascension = %g°" % right_ascension
        print "true_anomaly    = %g°" % true_anomaly
        print "_______________________________________________________________________________________"
        print "Current position of the satellite"
        print "Latitude        =",lat
        print "Longitude       =",lon
        print "----------------------------------------------------------------------------------------"
# Following is an example of two line element dateset format
# ISS (ZARYA)
# 1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
# 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537

names = ["stations"]

for i in names:
    f = urllib.urlopen("http://www.celestrak.com/NORAD/elements/%s.txt" % i)
    ele = ""
    for j in f:
        ele += j
        if (j[0] == '2'):
            ele = ele.strip()
            if ele.startswith("HAVELSAT"):
                print_this(ele, Info=True, labels=True)
            ele = ""
            



# End of File
