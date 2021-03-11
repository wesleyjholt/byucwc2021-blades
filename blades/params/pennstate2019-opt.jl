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
precone = 0.0*pi/180

# define hub, tip, and radial sections
Rhub = .0225
Rtip = .225
nsections = 17  # number of radial sections
r_spacing = (Rtip - Rhub)/nsections # distance between each radial section
r = collect(range(0, stop=1.0, length=nsections)).*(Rtip - Rhub - r_spacing) .+ Rhub .+ r_spacing/2 # evenly space the radial sections along the length of the blade
# r = [0.025570089, 0.03889375, 0.052217411, 0.064974107, 0.07858125, 0.091904911, 0.104945089, 0.11826875, 0.131592411, 0.144916071, 0.158239732, 0.171563393, 0.184603571, 0.197927232, 0.211250893, 0.224574554]

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

# chord and twist
chord = [0.1996870833343525, 0.19981302263435383, 0.1999167100825758, 0.19998331111116546, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
theta = [0.9705563982926866, 0.9705563982926866, 0.9648984007349249, 0.9590290663380072, 0.9565710572545105, 0.9566103305409216, 0.9582328432537264, 0.9636205047113439, 0.9786063905633794, 0.9982953741466316, 1.0171543993067278, 1.0296504098892958, 1.0302121657363474, 1.016444274968446, 0.9903126150073802, 0.9540059491279926, 0.9097130406051254]

# chord and twist spline parameters
cspline = [0.1996204050423886, 0.2, 0.2, 0.2, 0.2]
tspline = [0.9705563982926866, 0.9588916793356592, 1.0316277654620356, 0.8852557338757749]

# optimal tip-speed ratio
tsr = 4.0

# optimal pitch angles
pitch = [0.0, 0.0029589378125453036, 0.01148068785166434, 0.01525362800745561, 0.018470283855667124, 0.020323053224882093, 0.023174851728522794]