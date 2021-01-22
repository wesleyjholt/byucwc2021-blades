# This file sets the parameters for Penn State's 2019 blade for the CWC.

# import packages
using CCBlade
using DelimitedFiles

# import the radial section properties
sectiondata = readdlm("input-files/pennstate2019.txt", skipstart=1)

# separate radial position, chord, and twist
r = sectiondata[:,1]
chord = sectiondata[:,2]
theta = sectiondata[:,3]*pi/180

# Define airfoil(s). This blade only has one airfoil that is used throughout the entire length.
naf = 1
aftypes = Array{AlphaAF}(undef, naf)
aftypes[1] = AlphaAF("airfoil-data/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = Int.(sectiondata[:,4])

# create airfoil array
airfoils = aftypes[af_idx]

# define sections
sections = Section.(r, chord, theta, airfoils)

# define rotor
Rhub = r[1]
Rtip = r[end]
B = 3
precone = 0.0*pi/180
du = DuSeligEggers()  # Du-Selig correction for lift, Eggers correction for drag

rotor = Rotor(Rhub, Rtip, B, rotation=du, precone=precone, turbine=true)