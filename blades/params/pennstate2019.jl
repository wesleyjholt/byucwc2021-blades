# This file sets the parameters for Penn State's 2019 blade for the CWC.

# import packages
using CCBlade
using DelimitedFiles
using PyPlot
import FLOWMath

# import the radial section properties
sectiondata = readdlm("input-files/pennstate2019.txt", skipstart=1)

# separate radial position, chord, and twist
r = sectiondata[:,1]
chord = sectiondata[:,2]
theta = sectiondata[:,3]*pi/180
# figure()
# plot(chord, r)
# xlabel("radial position")
# ylabel("chord")
# plt.show()

# import Wortmann FX 63-137 airfoil data
# Re_vec = [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 
#     0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.13, 0.16, 0.2, 0.25, 0.3] .* 1e6
Re_vec = [0.005, 0.008, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.13, 0.16, 0.2, 0.25, 0.3] .* 1e6
airfoildata_directory = "airfoil-data/wortmann-fx63137/"
wortmann_fx_63137_filenames = Vector{String}(undef, length(Re_vec))
for i = 1:length(Re_vec)
    wortmann_fx_63137_filenames[i] = airfoildata_directory * "wortmann-fx63137-extrapolated-Re" * rpad(Re_vec[i]/1e6, 5, '0') * "e6.dat"
    alpha_vec = readdlm(wortmann_fx_63137_filenames[i], skipstart=3)[:,1]
end
wortmann_fx_63137 = AlphaReAF(wortmann_fx_63137_filenames; radians=true)

# Define airfoil(s). This blade only has one airfoil that is used throughout the entire length.
naf = length(r)
aftypes = Array{AlphaReAF}(undef, naf)
aftypes[1] = wortmann_fx_63137

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
