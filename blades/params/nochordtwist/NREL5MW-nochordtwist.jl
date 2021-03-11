#= 
DESCRIPTION: Defines rotor and blade parameters/files for the NREL 5MW wind turbine blades.
NOTES: This file is intended to be used by a blade analysis and/or optimization script.
=#

# define rotor parameters
B = 3   # number of blades
turbine = true
precone = 2.5*pi/180

# define hub, tip, and radial sections
Rhub = 1.5
Rtip = 63.0
r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
    28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
    56.1667, 58.9000, 61.6333]

# In this case we have 8 different airfoils that we load into an array. These airfoils are defined in files.
aftypes = Array{AlphaAF}(undef, 8)
aftypes[1] = AlphaAF("airfoil-data/example-airfoil-data/Cylinder1.dat", radians=false)
aftypes[2] = AlphaAF("airfoil-data/example-airfoil-data/Cylinder2.dat", radians=false)
aftypes[3] = AlphaAF("airfoil-data/example-airfoil-data/DU40_A17.dat", radians=false)
aftypes[4] = AlphaAF("airfoil-data/example-airfoil-data/DU35_A17.dat", radians=false)
aftypes[5] = AlphaAF("airfoil-data/example-airfoil-data/DU30_A17.dat", radians=false)
aftypes[6] = AlphaAF("airfoil-data/example-airfoil-data/DU25_A17.dat", radians=false)
aftypes[7] = AlphaAF("airfoil-data/example-airfoil-data/DU21_A17.dat", radians=false)
aftypes[8] = AlphaAF("airfoil-data/example-airfoil-data/NACA64_A17.dat", radians=false)

# aftypes = Array{AlphaReAF}(undef, 8)
# aftypes[1] = AlphaReAF(["airfoil-data/example-airfoil-data/Cylinder1.dat", "airfoil-data/example-airfoil-data/Cylinder1.dat"], radians=false)
# aftypes[2] = AlphaReAF(["airfoil-data/example-airfoil-data/Cylinder2.dat", "airfoil-data/example-airfoil-data/Cylinder2.dat"], radians=false)
# aftypes[3] = AlphaReAF(["airfoil-data/example-airfoil-data/DU40_A17.dat", "airfoil-data/example-airfoil-data/DU40_A17.dat"], radians=false)
# aftypes[4] = AlphaReAF(["airfoil-data/example-airfoil-data/DU35_A17.dat", "airfoil-data/example-airfoil-data/DU35_A17.dat"], radians=false)
# aftypes[5] = AlphaReAF(["airfoil-data/example-airfoil-data/DU30_A17.dat", "airfoil-data/example-airfoil-data/DU30_A17.dat"], radians=false)
# aftypes[6] = AlphaReAF(["airfoil-data/example-airfoil-data/DU25_A17.dat", "airfoil-data/example-airfoil-data/DU25_A17.dat"], radians=false)
# aftypes[7] = AlphaReAF(["airfoil-data/example-airfoil-data/DU21_A17.dat", "airfoil-data/example-airfoil-data/DU21_A17.dat"], radians=false)
# aftypes[8] = AlphaReAF(["airfoil-data/example-airfoil-data/NACA64_A17.dat", "airfoil-data/example-airfoil-data/NACA64_A17.dat"], radians=false)

# indices correspond to which airfoil is used at which station
af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]

# create airfoil array
airfoils = aftypes[af_idx]
