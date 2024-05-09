#geometric things

import numpy as np

from bibliotecas.argon_param import *

def solid_angle_on_axis(side1,side2,dist):
    
    return 4*np.arctan((side1*side2/(4*dist**2))/np.sqrt(1+(side1/(2*dist))**2+(side2/(2*dist))**2))


def solid_angle(x0,y0,z0, l1, l2, dev_x1, dev_y1, dev_z1): # divided by 4pi
    
    A= dev_y1-y0-l1
    
    B= dev_z1-z0-l2
    
    d = np.abs(x0-dev_x1)
    
    angle = (solid_angle_on_axis(2*np.abs(A+2*l1),2*np.abs(B+2*l2),d) - np.sign(A)*solid_angle_on_axis(2*np.abs(A),2*np.abs(B+2*l2),d) - np.sign(B)*solid_angle_on_axis(2*np.abs(A+2*l1),2*np.abs(B),d) + np.sign(A)*np.sign(B)*solid_angle_on_axis(2*np.abs(A),2*np.abs(B),d))/(16*np.pi)

    return angle


def distance(x, y, z, x0, y0, z0):
    return np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)


def semi_major(t):
    return v*t/2

def semi_minor_sqr(t, x, y, z, x0, y0, z0):
    value = (v**2*t**2 - distance( x, y, z, x0, y0, z0)**2)/4
    if value <0:
        value =0
    return value

def eccentricity(t, x, y, z, x0, y0, z0):
    return distance(x, y, z, x0, y0, z0)/(v*t)

def radial_param(a, e, phi):
    return (a*(1-e**2))/(1-e*np.cos(phi))

def vol_ellipsoid(t, x, y, z, x0, y0, z0):
    distance=np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
    if semi_minor_sqr(t, x, y, z, x0, y0, z0)< 0:
        value = 1
    else:
        value = 4/3*np.pi*semi_major(t)*semi_minor_sqr(t, x, y, z, x0, y0, z0)
    return value

def normal_cos(n1, n2, n3, v1, v2, v3):
    value = (n1*v1+n2*v2+n3*v3)/np.sqrt(v1**2+v2**2+v3**2)
    return value

def normal_sin(n1, n2, n3, v1, v2, v3):
    value = np.sqrt(1 - (n1*v1+n2*v2+n3*v3)**2/np.sqrt(v1**2+v2**2+v3**2)**2)
    return value


def cos_direction_angle(t, x, y, z, x0, y0, z0, n1, n2, n3,semi_minor_sqr,semi_major,normal_sin,normal_cos):
    vector = np.array([semi_minor_sqr*normal_sin/(semi_minor_sqr*normal_sin**2 + semi_major**2*normal_cos**2),
                       0,
                      1 - semi_minor_sqr*normal_sin**2/(normal_cos*(semi_minor_sqr*normal_sin**2 + semi_major**2*normal_cos))])
    value = normal_cos/np.sqrt(vector.dot(vector))
    return value

def int_cap_volume(t, x, y, z, x0, y0, z0, n1, n2, n3,semi_minor_sqr,semi_major,normal_sin,normal_cos):
    distance=np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
    if t*v - distance< 0:
        value = 0
    else: 
        value = np.pi*semi_minor_sqr*semi_major*(2/3 + distance*normal_cos/2/(semi_minor_sqr*normal_sin**2+semi_major**2*normal_cos**2)**0.5 - distance**3*normal_cos**3/24/(semi_minor_sqr*normal_sin**2+semi_major**2*normal_cos**2)**1.5)
    return value

def ext_cap_volume(t, x, y, z, x0, y0, z0, n1, n2, n3,semi_minor_sqr,semi_major,normal_sin,normal_cos):
    distance=np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
    if t*v - distance< 0:
        value = 0
    else:    
        value = np.pi*semi_minor_sqr*semi_major*(2/3 - distance*normal_cos/2/(semi_minor_sqr*normal_sin**2+semi_major**2*normal_cos**2)**0.5 + distance**3*normal_cos**3/24/(semi_minor_sqr*normal_sin**2+semi_major**2*normal_cos**2)**1.5)
    return value


