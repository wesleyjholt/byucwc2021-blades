# This file sets the parameters for Penn State's 2019 blade for the CWC.

# import packages
using CCBlade
using DelimitedFiles
using PyPlot
import FLOWMath

# define scale
scale = 1e0

# define the radial section properties
r = scale .* [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
    28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
    56.1667, 58.9000, 61.6333]
chord = scale .* [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
    3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
theta = pi/180 .* [13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
    6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

# Define airfoils.  In this case we have 8 different airfoils that we load into an array.
# These airfoils are defined in files.
aftypes = Array{AlphaAF}(undef, 8)
aftypes[1] = AlphaAF("airfoil-data/Cylinder1.dat", radians=false)
aftypes[2] = AlphaAF("airfoil-data/Cylinder2.dat", radians=false)
aftypes[3] = AlphaAF("airfoil-data/DU40_A17.dat", radians=false)
aftypes[4] = AlphaAF("airfoil-data/DU35_A17.dat", radians=false)
aftypes[5] = AlphaAF("airfoil-data/DU30_A17.dat", radians=false)
aftypes[6] = AlphaAF("airfoil-data/DU25_A17.dat", radians=false)
aftypes[7] = AlphaAF("airfoil-data/DU21_A17.dat", radians=false)
aftypes[8] = AlphaAF("airfoil-data/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]

# create airfoil array
airfoils = aftypes[af_idx]

# define sections
sections = Section.(r, chord, theta, airfoils)

# define rotor
Rhub = scale * 1.5
Rtip = scale * 63.0
B = 3
precone = 2.5*pi/180

rotor = Rotor(Rhub, Rtip, B, precone=precone, turbine=true)
