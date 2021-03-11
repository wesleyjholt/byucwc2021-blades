# This file sets the parameters for Penn State's 2019 blade for the CWC.

# import packages
using CCBlade
using DelimitedFiles
using PyPlot
import FLOWMath

# import the radial section properties
sectiondata = readdlm("input-files/testblades01.txt", skipstart=1)

# set scale
scale = 1e0

# separate radial position, chord, and twist
r = sectiondata[:,1]
chord = sectiondata[:,2]
theta = sectiondata[:,3]*pi/180

# r = scale .* [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
#     28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
#     56.1667, 58.9000, 61.6333]
# chord = scale .* [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
#     3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
# theta = pi/180 .* [13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
#     6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

# r = r_root_to_tip[2:end-1]
# chord = chord_root_to_tip[2:end-1]
# theta = theta_root_to_tip[2:end-1]
# figure()
# plot(chord, r)
# xlabel("radial position")
# ylabel("chord")
# plt.show()

# import USNPS4 airfoil data

# set airfoil data directory
af_directory = "airfoil-data/"

# specify name of airfoil
af_name = "NACA 4412"    # this should match exactly with the name of the directory that contains all the data for the airfoil of interest

# set the Reynolds numbers
Re_vec = [.005, .01, .06, .08, .1, .13, .15] .* 1e6    # you should already have extrapolated airfoil data for these Re numbers, saved in the appropriate location

# save 2D extrapolated airfoil data file
af_files = Vector{String}(undef, length(Re_vec))
for i = 1:length(Re_vec)
    af_files[i] = af_name * "/extrapolated/" * af_name * "_Re" * rpad(Re_vec[i]/1e6, 5, '0') * "e6.dat"
end
af = AlphaReAF(af_files; radians=true)

airfoildata_directory = "airfoil-data/usnps4/smooth-extrapolated/"
usnps4_filenames = Vector{String}(undef, length(Re_vec))
for i = 1:length(Re_vec)
    usnps4_filenames[i] = airfoildata_directory * "USNPS4-extrapolated-Re" * rpad(Re_vec[i]/1e6, 5, '0') * "e6.dat"
    alpha_vec = readdlm(usnps4_filenames[i], skipstart=3)[:,1]
end
usnps4 = AlphaReAF(usnps4_filenames; radians=true)
usnps4_oneRe = AlphaAF(airfoildata_directory * "USNPS4-extrapolated-Re" * rpad(Re_vec[1]/1e6, 5, '0') * "e6.dat", radians=true)

# import example airfoil data
# du30_a17 = AlphaAF("airfoil-data/DU30_A17.dat", radians=false)

# # import USNPS4 airfoil data (extrapolated by QBlade)
# Re_vec = [.005, .01, .02, .03, .04, .05, .06, .08, .1, .13] .* 1e6
# airfoildata_directory = "airfoil-data/usnps4-extrapolated-qblade/"
# usnps4_filename = airfoildata_directory * "USNPS4_T1_Re" * rpad(Re_vec[1]/1e6, 5, '0') * "_M0.00_N9.0_360_V.dat"
# alpha_vec = readdlm(usnps4_filename, skipstart=14)[:,1] .* pi/180
# cl = zeros(length(alpha_vec), length(Re_vec))
# cd = zeros(length(alpha_vec), length(Re_vec))
# usnps4_filenames = Vector{String}(undef, length(Re_vec))
# for i = 1:length(Re_vec)
#     usnps4_filename = airfoildata_directory * "USNPS4_T1_Re" * rpad(Re_vec[1]/1e6, 5, '0') * "_M0.00_N9.0_360_V.dat"
#     data = readdlm(usnps4_filename, skipstart=14)
#     cl[:,i] = data[:,2]
#     cd[:,i] = data[:,3]
# end
# usnps4 = AlphaReAF(alpha_vec, Re_vec, cl, cd)

# Define airfoil(s). This blade only has one airfoil that is used throughout the entire length.
naf = length(r)
aftypes = Array{AlphaReAF}(undef, naf)
# aftypes = Array{AlphaAF}(undef, naf)
aftypes[1] = usnps4
# aftypes[1] = usnps4_oneRe

# indices correspond to which airfoil is used at which station
af_idx = Int.(sectiondata[:,4])
# af_idx = ones(Int,length(r))

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
