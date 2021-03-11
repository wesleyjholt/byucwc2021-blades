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
chord = [3.6619356158414518, 3.8634044265178944, 4.0065197560086485, 4.0902609096941225, 4.029026578309546, 3.7481440477515515, 3.2826733610785257, 2.8010678999524115, 2.475636603945036, 2.3535411838953335, 2.3041038694666898, 2.232874378467609, 2.050156756743517, 1.7888312310867889, 1.5322873310573581, 1.298702548082692, 1.037620410270333]
theta = [0.2776153699751321, 0.2776153699751321, 0.2776148105152295, 0.22417513309782494, 0.16577062883613677, 0.11802779320703864, 0.08562554080851306, 0.07104711259249598, 0.06303957899585656, 0.05703846402039875, 0.05036854283788042, 0.040354590620059375, 0.02798092485069706, 0.01562902530882623, 0.00493422814846615, -0.004125627857854169, -0.01381966272315736]

# chord and twist spline parameters
cspline = [3.542, 3.9830040205603026, 2.475636603945036, 2.1083333964631272, 0.896100804365593]
tspline = [0.2776153699751321, 0.07663929413258447, 0.040354590620059375, -0.018958342891837683]

# optimal tip-speed ratio
tsr = 9.103069195672084

# optimal pitch angles
pitch = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]