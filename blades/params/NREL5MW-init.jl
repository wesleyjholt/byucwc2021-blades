# ================================================
# ============ MACHINE GENERATED FILE ============
# ================================================


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


# chord and twist
chord = [3.7204206406588893, 4.059114731182795, 4.348522281903936, 4.595075940860216, 4.644277768817204, 4.442398193843963, 4.219073181479043, 3.9871806990546057, 3.748, 3.50700594073793, 3.2653808915282188, 3.017965376555897, 2.7574272038717482, 2.4375773341399474, 2.1239594090733505, 1.8509418882839281, 1.564943233910539]
theta = [0.23226841685540536, 0.23226841685540536, 0.23226820766473041, 0.21201226319794977, 0.18860594202739364, 0.16696488790544434, 0.14791797293307937, 0.13187925797072722, 0.11675964214737788, 0.10187658330568483, 0.08691467213677301, 0.0715584993317675, 0.0559307715223736, 0.04031679117583097, 0.027263919805328163, 0.01676515687104599, 0.006193043980894238]

# chord and twist spline parameters
cspline = [3.542, 4.6035, 3.748, 2.8255, 1.419]
tspline = [0.23226841685540536, 0.13962634015954636, 0.0715584993317675, 0.0008726646259971648]

# optimal tip-speed ratio
tsr = 7.55

# optimal pitch angles
pitch = [0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659, 0.3490658503988659]