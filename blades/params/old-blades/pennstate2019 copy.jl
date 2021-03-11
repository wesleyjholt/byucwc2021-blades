#= 
DESCRIPTION: Defines rotor and blade parameters/files for the NREL 5MW wind turbine blades.
NOTES: This file is intended to be used by a blade analysis and/or optimization script.
=#

# define rotor parameters
B = 3   # number of blades
turbine = true
precone = 0.0*pi/180

# define hub, tip, and radial sections
Rhub = .0225
Rtip = .225
r = [0.025570089, 0.03889375, 0.052217411, 0.064974107, 0.07858125, 0.091904911, 0.104945089, 0.11826875, 0.131592411, 0.144916071, 0.158239732, 0.171563393, 0.184603571, 0.197927232, 0.211250893, 0.224574554]

# get the airfoil data
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

# create airfoil array
af_idx = Int.(ones(length(r)))
airfoils = aftypes[af_idx]