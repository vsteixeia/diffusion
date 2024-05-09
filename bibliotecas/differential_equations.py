#differential equations solutions

import numpy as np

from bibliotecas.argon_param import *
import bibliotecas.geofunc2 as g2
from scipy import integrate

from scipy import special


#from Fourier-Transf-Desperado-3D-Try8 import semi_major, semi_minor_sqr, eccentricity, vol_ellipsoid, normal_cos, normal_sin, cos_direction_angle, int_cap_volume, ext_cap_volume

def distance(x, y, z, x0, y0, z0):
    return np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)


def delta_fun(x, a=5):
    return 1/a*(np.heaviside(x,1) - np.heaviside(x-a,1))



def solution_wave(t, x, y, z, x0, y0, z0):
    value = np.exp(-(gamma+alpha)*t)*((8-3*np.exp(-gamma*t)+2*gamma*t+4*gamma**2*t**2)
                    *delta_fun(v*t-distance(x, y, z, x0, y0, z0))/distance(x, y, z, x0, y0, z0)**2
            )
    return S/(20*np.pi)*value

#"""
def deposit_wave_point(x, y, z, x0, y0, z0): #multiplied by v
    value = S/(20*np.pi)*(
        np.exp(-(gamma+alpha)*(distance(x, y, z, x0, y0, z0)/v))*((8
    -3*np.exp(-gamma*(distance(x, y, z, x0, y0, z0)/v))+2*gamma*(distance(x, y, z, x0, y0, z0)/v)
    +4*gamma**2*(distance(x, y, z, x0, y0, z0)/v)**2)/distance(x, y, z, x0, y0, z0)**2
                                                                 ))*(x0-x)/np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
    return value


def deposit_wave(x, y, z, x0, y0, z0): #multiplied by v
    aux = lambda eps2, eps3: deposit_wave_point(x, y+eps2, z+eps3, x0, y0, z0)                                                   
    value = integrate.dblquad(aux, -l, l, -l, l)[0]
    return value
#"""

def solution_dif(t, x, y, z, x0, y0, z0):
    if v*t >= distance(x, y, z, x0, y0, z0):
        value = S/(20*np.pi)*np.exp(-(gamma+alpha)*t)*(gamma**2/v*(1/(v*np.sqrt(v**2*t**2-distance(x, y, z, x0, y0, z0)**2))*special.iv(1, gamma*np.sqrt(v**2*t**2-distance(x, y, z, x0, y0, z0)**2)/v)
                                                 +4*t/(v**2*t**2-distance(x, y, z, x0, y0, z0)**2)*special.iv(2, gamma*np.sqrt(v**2*t**2-distance(x, y, z, x0, y0, z0)**2)/v)
                                                )
                )
    else:
        value = 0
    return value

def solution(t, x, y, z, x0, y0, z0):
    return solution_wave(t, x, y, z, x0, y0, z0) + solution_dif(t, x, y, z, x0, y0, z0)



def excess(t, x, y, z, x0, y0, z0):

    thing = 0
    #backscattering on x
    if  t-distance(x, y, z, -x0, y0, z0)/v > 0:
        x00 = 0
        y00 = x0/(x+x0)*(y-y0)+y0
        z00 = x0/(x+x0)*(z-z0)+z0
        t0 = distance(x00, y00, z00, x, y, z)/v
        semi_major0 = g2.semi_major(t-t0)
        semi_minor_sqr0 = g2.semi_minor_sqr(t-t0, x00, y00, z00, x0, y0, z0)
        vol_ellipsoid0 = g2.vol_ellipsoid(t-t0, x00, y00, z00, x0, y0, z0)
        normal_cos0 = g2.normal_cos(1, 0, 0, x0-x00, y0-y00, z0-z00)
        normal_sin0 = g2.normal_sin(1, 0, 0, x0-x00, y0-y00, z0-z00)
        ext_cap_volume0 = g2.ext_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 1, 0, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        int_cap_volume0 = g2.int_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 1, 0, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        thing += ext_cap_volume0/int_cap_volume0*solution_dif(t, x, y, z, -x0, y0, z0)
    if t - distance(x, y, z, 2*Lx-x0, y0, z0)/v > 0:
        x00 = Lx
        y00 = (Lx-x0)/(2*Lx -x-x0)*(y-y0)+y0
        z00 = (Lx-x0)/(2*Lx -x-x0)*(z-z0)+z0
        t0 = distance(x00, y00, z00, x, y, z)/v
        semi_major0 = g2.semi_major(t-t0)
        semi_minor_sqr0 = g2.semi_minor_sqr(t-t0, x00, y00, z00, x0, y0, z0)
        vol_ellipsoid0 = g2.vol_ellipsoid(t-t0, x00, y00, z00, x0, y0, z0)
        normal_cos0 = g2.normal_cos(-1, 0, 0, x0-x00, y0-y00, z0-z00)
        normal_sin0 = g2.normal_sin(-1, 0, 0, x0-x00, y0-y00, z0-z00)
        ext_cap_volume0 = g2.ext_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, -1, 0, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        int_cap_volume0 = g2.int_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, -1, 0, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        thing += ext_cap_volume0/int_cap_volume0*solution_dif(t, x, y, z, 2*Lx-x0, y0, z0)
    #backscateriing on y
    if  t-distance(x, y, z, x0, -y0, z0)/v > 0:
        x00 = y0/(y+y0)*(x-x0)+x0
        y00 = 0
        z00 = y0/(y+y0)*(z-z0)+z0
        t0 = distance(x00, y00, z00, x, y, z)/v
        semi_major0 = g2.semi_major(t-t0)
        semi_minor_sqr0 = g2.semi_minor_sqr(t-t0, x00, y00, z00, x0, y0, z0)
        vol_ellipsoid0 = g2.vol_ellipsoid(t-t0, x00, y00, z00, x0, y0, z0)
        normal_cos0 = g2.normal_cos(0, 1, 0, x0-x00, y0-y00, z0-z00)
        normal_sin0 = g2.normal_sin(0, 1, 0, x0-x00, y0-y00, z0-z00)
        ext_cap_volume0 = g2.ext_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, 1, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        int_cap_volume0 = g2.int_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, 1, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        thing += ext_cap_volume0/int_cap_volume0*solution_dif(t, x, y, z, x0, -y0, z0)
    if t - distance(x, y, z, x0, 2*Ly-y0, z0)/v > 0:
        x00 = (Ly-y0)/(2*Ly -y-y0)*(x-x0)+x0
        y00 = Ly
        z00 = (Ly-y0)/(2*Ly -y-y0)*(z-z0)+z0
        t0 = distance(x00, y00, z00, x, y, z)/v
        semi_major0 = g2.semi_major(t-t0)
        semi_minor_sqr0 = g2.semi_minor_sqr(t-t0, x00, y00, z00, x0, y0, z0)
        vol_ellipsoid0 = g2.vol_ellipsoid(t-t0, x00, y00, z00, x0, y0, z0)
        normal_cos0 = g2.normal_cos(0, -1, 0, x0-x00, y0-y00, z0-z00)
        normal_sin0 = g2.normal_sin(0, -1, 0, x0-x00, y0-y00, z0-z00)
        ext_cap_volume0 = g2.ext_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, -1, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        int_cap_volume0 = g2.int_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, -1, 0,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        thing += ext_cap_volume0/int_cap_volume0*solution_dif(t, x, y, z, x0, 2*Ly-y0, z0)
    #backscateriing on z
    if  t-distance(x, y, z, x0, y0, -z0)/v > 0:
        x00 = z0/(z+z0)*(x-x0)+x0
        y00 = z0/(z+z0)*(y-y0)+y0
        z00 = 0
        t0 = distance(x00, y00, z00, x, y, z)/v
        semi_major0 = g2.semi_major(t-t0)
        semi_minor_sqr0 = g2.semi_minor_sqr(t-t0, x00, y00, z00, x0, y0, z0)
        vol_ellipsoid0 = g2.vol_ellipsoid(t-t0, x00, y00, z00, x0, y0, z0)
        normal_cos0 = g2.normal_cos(0, 0, 1, x0-x00, y0-y00, z0-z00)
        normal_sin0 = g2.normal_sin(0, 0, 1, x0-x00, y0-y00, z0-z00)
        ext_cap_volume0 = g2.ext_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, 0, 1,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        int_cap_volume0 = g2.int_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, 0, 1,semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        thing += ext_cap_volume0/int_cap_volume0*solution_dif(t, x, y, z, x0, y0, -z0)
    if t - distance(x, y, z, x0, y0, 2*Lz-z0)/v > 0:
        x00 = (Lz-z0)/(2*Lz -z-z0)*(x-x0)+x0
        y00 = (Lz-z0)/(2*Lz -z-z0)*(y-y0)+y0
        z00 = Lz
        t0 = distance(x00, y00, z00, x, y, z)/v
        semi_major0 = g2.semi_major(t-t0)
        semi_minor_sqr0 = g2.semi_minor_sqr(t-t0, x00, y00, z00, x0, y0, z0)
        vol_ellipsoid0 = g2.vol_ellipsoid(t-t0, x00, y00, z00, x0, y0, z0)
        normal_cos0 = g2.normal_cos(0, 0, -1, x0-x00, y0-y00, z0-z00)
        normal_sin0 = g2.normal_sin(0, 0, -1, x0-x00, y0-y00, z0-z00)
        ext_cap_volume0 = g2.ext_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, 0, -1, semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        int_cap_volume0 = g2.int_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, 0, 0, -1, semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        thing += ext_cap_volume0/int_cap_volume0*solution_dif(t, x, y, z, x0, y0, 2*Lz-z0)

    return -thing




def draw2d(t, x, y, z ,x0,y0,z0):
    a = len(x)
    b = len(y)
    value = np.zeros((a,b))
    for i in range(a):
        for j in range(b):
            value[i,b-j-1] = S/(20*np.pi)*(
                excess(t, x[i], y[j], z, x0, y0, z0)
                #+solution_wave(t, x[i], y[j], z0, x0, y0, z0)
                +solution_dif(t, x[i], y[j], z, x0, y0, z0)
            )
    #print(np.max(-value))
    return np.transpose(value)




def flux(eps2, eps3, t, x, y, z, x0, y0, z0, n1, n2, n3):
    if v*t > distance(x, eps2 +y, eps3+z, x0, y0, z0):
        x00 = 0
        y00 = x0/(x+x0)*(y-y0)+y0
        z00 = x0/(x+x0)*(z-z0)+z0
        t0 = distance(x00, y00, z00, x, y, z)/v
        semi_major0 = g2.semi_major(t-t0)
        semi_minor_sqr0 = g2.semi_minor_sqr(t-t0, x00, y00, z00, x0, y0, z0)
        vol_ellipsoid0 = g2.vol_ellipsoid(t-t0, x00, y00, z00, x0, y0, z0)
        normal_cos0 = g2.normal_cos(n1, n2, n3, x0-x00, y0-y00, z0-z00)
        normal_sin0 = g2.normal_sin(n1, n2, n3, x0-x00, y0-y00, z0-z00)
        ext_cap_volume0 = g2.ext_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, n1, n2, n3, semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        int_cap_volume0 = g2.int_cap_volume(t-t0, x00, y00, z00, x0, y0, z0, n1, n2, n3, semi_minor_sqr0,semi_major0,normal_sin0,normal_cos0)
        

        thing = int_cap_volume0/vol_ellipsoid0*normal_cos0
        value = v*thing*(
                excess(t, x, y+eps2, z+eps3, x0, y0, z0)
                +solution_wave(t, x, y+eps2, z+eps3, x0, y0, z0)
                +solution_dif(t, x, y+eps2, z+eps3, x0, y0, z0)
            )
    else:
        value = 0
    return value

"""
def flux_dif(eps2, eps3, t, n1, n2, n3):
    if v*t > distance(dev_x1, eps2 + dev_y1, eps3 + dev_z1, x0, y0, z0):
        value = v*(average_outward3(t, dev_x1, eps2 + dev_y1, eps3 + dev_z1, x0, y0, z0, n1, n2, n3)*(
                +(solution_dif(t, dev_x1, eps2 + dev_y1, eps3 + dev_z1, x0, y0, z0)*(x0-dev_x1)/np.sqrt((dev_x1-x0)**2+(dev_y1+eps2-y0)**2+(dev_z1+eps3-z0)**2)
                -excess(t, dev_x1, eps2 + dev_y1, eps3 + dev_z1, x0, y0, z0)))
                )
    else:
        value = 0
    return value
"""

def flux_detector(t_max, x, y, z, x0, y0, z0, n1, n2, n3, skip = 0.1):
    abs_tol = 1
    dim = 2
    
    t = np.arange(0, t_max, skip)
    value = []
    
    N = 500
    
    for i in range(len(t)):
        accum = 0
        for j in range(N):
            eps2 = np.random.uniform(-l, l)
            eps3 = np.random.uniform(-l, l)
            accum += flux(eps2, eps3, t[i], x, y, z, x0, y0, z0, n1, n2, n3)
        value.append(accum)
    measure = 4 * l**2/N
    value = measure * np.array(value)
    
    return value


def total_flow_detector(t_max, x, y, z, x0, y0, z0, n1, n2, n3, skip = 0.1):
    
    spread_data = flux_detector(t_max, x, y, z, x0, y0, z0, n1, n2, n3, skip = skip)
    value = sum(spread_data)*skip*t_max
    
    return value

"""
def Total_flow(t_max, dev_x1, dev_y1, dev_z1, x0, y0, z0, n1, n2, n3 skip = 0.1):
    dim = 2
    
    t = np.linspace(0, t_max, num_tikz)
    value = []
    
    N = 5_000
    accum = 0
    for i in range(N):
        t = np.random.uniform(0, t_max)
        eps2 = np.random.uniform(-l, l)
        eps3 = np.random.uniform(-l, l)
        accum += flux_dif(eps2, eps3, t, n1, n2, n3)
    value.append(accum)
    measure = 4 * l**2/N
    value = measure * np.array(value)
    
    return value
"""
