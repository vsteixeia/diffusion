import numpy as np

# constants # Andrzej # Galimov

lambda_abs = 2000 # absorption length
lambda_abs1 = 2000

lambda_rs =  100#100 # scattering lenght

g = 0.025 # average scattering cossine #anisotropy for Rayleigh scattering

lambda_rs_red = lambda_rs/(1-g) #100 # scattering lenght

v = 21.7 # 21.7 in one source cm/ns (should be phase velocity) # group speed of light in LAr # **13.4** A measurement of the group velocity of scintillation light in liquid argon
# using the index of refraction of the CERN paper we get a velocity of 22.08 ~  299792458÷1,358


Delta_E = 1

S = 100_000



# PD size # Andrzej

l = 9.3/2 # square



# Diffusion constants ????


alpha = -v/2*(2/(lambda_abs)+3/lambda_rs_red) + v/2*np.sqrt((2/(lambda_abs)+3/(lambda_rs_red))**2-(1/lambda_abs**2+3/(lambda_abs*lambda_rs_red)))

gamma = (alpha + v/2*(2/(lambda_abs)+3/lambda_rs_red))

# enlarging the detector

ed = 0


# Detector size # Andrzej

Lx = 1400+2*ed #cm

Ly = 365+2*ed #cm

Lz = 1200+2*ed #cm


# Flow through plane x = c

c = ed